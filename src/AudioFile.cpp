#include <iostream>

#include <vector>
#include <string>

extern "C" {
#include <libavformat/avformat.h>
};

#include "SampleSorter/Units.hpp"
#include "SampleSorter/AudioFile.hpp"
#include "SampleSorter/ProcessingException.hpp"

std::vector<std::vector<double> > AudioFile::read(
    std::string fileName,
    double startSeconds,
    double endSeconds,
    long * sampleRate) {

  if (startSeconds < 0 or endSeconds < 0) {
    throw ProcessingException("Time is negative!");
  }

  if (startSeconds > endSeconds) {
    throw ProcessingException("The time to start reading is after the time to end reading");
  }

  // Initialize FFmpeg
  av_register_all();

  // open file and get format information
  AVFormatContext * formatContext = NULL;
  if (avformat_open_input(&formatContext, fileName.c_str(), NULL, 0) != 0) {
    throw ProcessingException("Could not open audio file!");
  }

  // get stream info
  if (avformat_find_stream_info(formatContext, NULL) < 0) {
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not get stream info!");
  }

  // Find an audio stream
  // And its decoder
  AVCodec * codec = NULL;
  int audioStreamIndex = av_find_best_stream(formatContext, AVMEDIA_TYPE_AUDIO, -1, -1, &codec, 0);
  if (audioStreamIndex < 0) {
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not find audio stream or fetch its decoder!");
  }

  // Allocate context for decoding the codec
  AVCodecContext * codecContext = avcodec_alloc_context3(codec);
  if (codecContext == NULL) {
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not allocate a decoding context.");
  }

  // Fill the codecContext with parameters of the codec
  if (avcodec_parameters_to_context(codecContext, formatContext -> streams[audioStreamIndex] -> codecpar) != 0) {
    avcodec_close(codecContext);
    avcodec_free_context(&codecContext);
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not set codec context parameters.");
  }

  // Ask for non planar data
  codecContext -> request_sample_fmt = av_get_alt_sample_fmt(codecContext -> sample_fmt, 0);

  // Initialize the decoder
  if (avcodec_open2(codecContext, codec, NULL) != 0) {
    avcodec_close(codecContext);
    avcodec_free_context(&codecContext);
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not initialize the decoder");
  }

  // Allocate a frame
  AVFrame * frame = NULL;
  if ((frame = av_frame_alloc()) == NULL) {
    avcodec_close(codecContext);
    avcodec_free_context(&codecContext);
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not allocate audio frame!");
  }

  // prepare a packet
  AVPacket packet;
  av_init_packet(&packet);

  // END OF FFMPEG SET UP

  // fetch the sample rate
  *sampleRate = codecContext -> sample_rate;

  // Get start and end values in samples
  long startSample = Units::secondsToSamples(startSeconds, *sampleRate);
  long endSample = Units::secondsToSamples(endSeconds, *sampleRate);

  // Allocate the output vector
  std::vector<std::vector<double> > output(codecContext -> channels);
  for (long channel = 0; channel < codecContext -> channels; channel++) {
    output[channel].resize(endSample - startSample);
  }

  long sampleOffset = 0;
  int readError;

  // Read the file until either nothing is left
  // or we reach desired end of sample
  while (sampleOffset < endSample and
      (readError = av_read_frame(formatContext, &packet)) != AVERROR_EOF) {

    // Are we reading correctly?
    if (readError != 0) {
      throw ProcessingException("Read error!");
      break;
    }

    // Is this the correct stream?
    if (packet.stream_index != audioStreamIndex) {
      // Free the frame buffer and reset
      av_packet_unref(&packet);
      continue;
    }

    // Send the packet to the decoder!
    if (avcodec_send_packet(codecContext, &packet) == 0) {
      // Success!
      // Destroy it :)
      // free the buffers and reset fields
      av_packet_unref(&packet);
    } else  {
      // Failure!
      throw ProcessingException("Send packet eror!");
      break;
    }

    // receive the frame
    while ((readError = avcodec_receive_frame(codecContext, frame)) == 0) {
      // read the frame
      readFrame(codecContext, 
                frame, 
                output, 
                startSample, 
                endSample, 
                &sampleOffset);


      // free buffers and set fields to defaults
      av_frame_unref(frame);
    }

    if (readError != AVERROR(EAGAIN)) {
      throw ProcessingException("Receive packet eror!");
      break;
    }
  }

  // Drain the decoder

  // Clean up
  // Free all frame data
  av_frame_free(&frame);
  // Close the codec context
  avcodec_close(codecContext);
  // Free the context
  avcodec_free_context(&codecContext);
  avformat_close_input(&formatContext);

  return output;
}

void AudioFile::readFrame(
    const AVCodecContext * codecContext, 
    const AVFrame * frame,
    std::vector<std::vector<double> > & output,
    long startSample,
    long endSample,
    long * sampleOffset) {

  long numSamples = frame -> nb_samples;
  long numChannels = frame -> channels;

  // for every channel and sample in the frame
  for (long channel = 0; channel < numChannels; channel++) {
    for (long sample = 0; sample < numSamples; sample++) {

      // If the sample is within our sample range
      if (sample + *sampleOffset >= startSample
          and sample + *sampleOffset < endSample) {

        // get the sample value
        double value = getSample(codecContext, frame, sample, channel);

        // put it in the output
        output[channel][sample + *sampleOffset - startSample] = value;
      }
    }
  }

  *sampleOffset += numSamples;

}

double AudioFile::getSample(
    const AVCodecContext * codecContext,
    const AVFrame * frame,
    int sample,
    int channel) {

  // where it is stored
  uint8_t * frameBuffer;
  // its index
  int sampleLocation;

  // Is it planar?
  // I.E. Is each channel in its own vector
  if (av_sample_fmt_is_planar(codecContext -> sample_fmt)) {
    // yes, we search through the channel
    frameBuffer = frame -> extended_data[channel];
    sampleLocation = sample;
  } else {
    // no, we un-interleave them
    frameBuffer = frame -> extended_data[0];
    sampleLocation = sample * (codecContext -> channels) + channel;
  }

  // the number of bytes in an audio sample
  int sampleBytes = av_get_bytes_per_sample(codecContext -> sample_fmt);

  // the raw output value
  int64_t rawValue;

  switch(sampleBytes) {
    case 1:
      rawValue = ((uint8_t *)(frameBuffer))[sampleLocation];
      // its unsigned
      rawValue -= 127;
      break;

    case 2:
      rawValue = ((int16_t *)(frameBuffer))[sampleLocation];
      break;

    case 4:
      rawValue = ((int32_t *)(frameBuffer))[sampleLocation];
      break;

    case 8:
      rawValue = ((int64_t *)(frameBuffer))[sampleLocation];
      break;

    default:
      throw ProcessingException("Sample format is invalid!");
      return 0;
  }

  // if the raw values are padded for uneven bit depths (i.e. 24)
  rawValue = (rawValue >> (sampleBytes * 8 - codecContext -> bits_per_raw_sample));

  double value;

  switch (codecContext -> sample_fmt) {
    case AV_SAMPLE_FMT_U8:
    case AV_SAMPLE_FMT_U8P:
    case AV_SAMPLE_FMT_S16:
    case AV_SAMPLE_FMT_S16P:
    case AV_SAMPLE_FMT_S32:
    case AV_SAMPLE_FMT_S32P:
      // value / (2 ^ (number of Unsigned Bits) - 1)
      value = rawValue / ((double)( 1 << (sampleBytes * 8 - 1) - 1));
      break;
    case AV_SAMPLE_FMT_FLT:
    case AV_SAMPLE_FMT_FLTP:
      value = *((float *)(&rawValue));
      break;
    case AV_SAMPLE_FMT_DBL:
    case AV_SAMPLE_FMT_DBLP:
      value = *((double *)(&rawValue));
      break;
    default:
      throw ProcessingException("Sample format is invalid!");
      return 0;
  }

  return value;

}
