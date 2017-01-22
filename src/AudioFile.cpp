#include <iostream>

#include <vector>
#include <string>

extern "C" {
#include <libavformat/avformat.h>
};

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

  // Allocate a frame
  AVFrame * frame = av_frame_alloc();
  if (!frame) {
    throw ProcessingException("Could not allocate audio frame!");
  }

  // open file
  AVFormatContext * formatContext = NULL;
  if (avformat_open_input(&formatContext, fileName.c_str(), NULL, NULL) != 0) {
    av_free(frame);
    throw ProcessingException("Could not open audio file!");
  }

  // get stream info
  if (avformat_find_stream_info(formatContext, NULL) < 0) {
    av_free(frame);
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not get stream info!");
  }

  // Get the audio stream
  AVCodec * cdc = nullptr;
  int streamIndex = av_find_best_stream(formatContext, AVMEDIA_TYPE_AUDIO, -1, -1, &cdc, 0);
  if (streamIndex < 0) {
    // Could not find audio stream in file
    av_free(frame);
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not get stream info!");
  }

  AVStream * audioStream = formatContext -> streams[streamIndex];
  AVCodecContext * codecContext = audioStream -> codec;
  codecContext -> codec = cdc;

  if (avcodec_open2(codecContext, codecContext -> codec, NULL) != 0) {
    av_free(frame);
    avformat_close_input(&formatContext);
    throw ProcessingException("Could not open the context with decoder!");
  }

  *sampleRate = codecContext -> sample_rate;

  long startSample = startSeconds * (*sampleRate);
  long endSample = endSeconds * (*sampleRate);

  std::vector<std::vector<double> > output(codecContext -> channels);
  for (long channel = 0; channel < codecContext -> channels; channel++) {
    output[channel].resize(endSample - startSample);
  }

  long sampleOffset = 0;

  AVPacket readingPacket;
  av_init_packet(&readingPacket);

  while (sampleOffset < endSample and
      av_read_frame(formatContext, &readingPacket) == 0) {

    if (readingPacket.stream_index ==audioStream -> index) {
      AVPacket decodingPacket = readingPacket;

      while (decodingPacket.size > 0) {
        int gotFrame = 0;
        int result = avcodec_decode_audio4(codecContext, frame, &gotFrame, &decodingPacket);

        if (result >= 0 && gotFrame) {
          decodingPacket.size -= result;
          decodingPacket.data += result;

          if (sampleOffset + frame -> nb_samples > startSample) {
            readFrame(codecContext, 
                      frame, 
                      output, 
                      startSample, 
                      endSample, 
                      sampleOffset);
          }

          sampleOffset += frame -> nb_samples;

        } else {
          decodingPacket.size = 0;
          decodingPacket.data = nullptr;
        }
      }
    }
    av_packet_unref(&readingPacket);
  }

  // Clean up
  av_free(frame);
  avcodec_close(codecContext);
  avformat_close_input(&formatContext);

  return output;
}


void AudioFile::readFrame(
    const AVCodecContext * codecContext, 
    const AVFrame * frame,
    std::vector<std::vector<double> > & output,
    long startSample,
    long endSample,
    long sampleOffset) {

  long numSamples = frame -> nb_samples;
  long channels = frame -> channels;

  for (long channel = 0; channel < channels; channel++) {
    for (long sample = 0; sample < numSamples; sample++) {
      if (sample + sampleOffset >= startSample
          and sample + sampleOffset < endSample) {

        double value;

        switch(codecContext -> sample_fmt) {
          case AV_SAMPLE_FMT_U8:
            value = ((uint8_t *)(frame -> extended_data))[sample*channels + channel];
            break;
          case AV_SAMPLE_FMT_S16:
            value = ((int16_t *)(frame -> extended_data))[sample*channels + channel];
            break;
          case AV_SAMPLE_FMT_S32:
            value = ((int32_t *)(frame -> extended_data))[sample*channels + channel];
            break;
          case AV_SAMPLE_FMT_FLT:
            value = ((float *)(frame -> extended_data))[sample*channels + channel];
            break;
          case AV_SAMPLE_FMT_DBL:
            value = ((double *)(frame -> extended_data))[sample*channels + channel];
            break;
          case AV_SAMPLE_FMT_U8P:
            value = ((uint8_t *)(frame -> extended_data[channel]))[sample];
            break;
          case AV_SAMPLE_FMT_S16P:
            value = ((int16_t *)(frame -> extended_data[channel]))[sample];
            break;
          case AV_SAMPLE_FMT_S32P:
            value = ((int32_t *)(frame -> extended_data[channel]))[sample];
            break;
          case AV_SAMPLE_FMT_FLTP:
            value = ((float *)(frame -> extended_data[channel]))[sample];
            break;
          case AV_SAMPLE_FMT_DBLP:
            value = ((double *)(frame -> extended_data[channel]))[sample];
            break;
          default:
            throw ProcessingException("Sample format is invalid!");
            break;
        }

        output[channel][sample + sampleOffset - startSample] = value;
      }
    }
  }

}
