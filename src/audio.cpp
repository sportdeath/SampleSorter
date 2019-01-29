#include <vector>
#include <string>
#include <stdexcept>
#include <ciso646>

extern "C" {
#include <libavformat/avformat.h>
};

#include "sample_sorter/units.hpp"
#include "sample_sorter/audio.hpp"

using namespace sample_sorter;

void Audio::write(
    const std::vector<std::vector<double>> & audio,
    std::string filename) {

}

std::vector<std::vector<double>> Audio::read(
    std::string filename,
    double * sample_rate,
    double start_seconds,
    double end_seconds) {

  if (start_seconds < end_seconds) {
    throw std::invalid_argument(
        "The start read time is after the end read time!");
  }

  if (start_seconds < 0 or end_seconds < 0) {
    throw std::invalid_argument(
        "Start or end time is negative!");
  }

  // Begin to set up FFMPEG

  // Open the file and get format information
  AVFormatContext * format_context = NULL;
  if (avformat_open_input(&format_context, filename.c_str(), NULL, 0) != 0) {
    throw std::invalid_argument(
        "Could not open audio file: " + filename);
  }

  // Get stream info
  if (avformat_find_stream_info(format_context, NULL) < 0) {
    avformat_close_input(&format_context);
    throw std::runtime_error(
        "Could not get information about the stream in file: " + filename);
  }

  // Find an audio stream
  // And its decoder
  AVCodec * codec = NULL;
  int audio_stream_index = av_find_best_stream(
      format_context, 
      AVMEDIA_TYPE_AUDIO, 
      -1, -1, &codec, 0);
  if (audio_stream_index < 0) {
    avformat_close_input(&format_context);
    throw std::runtime_error(
        "Could not determine the best stream to use in the file: " + filename);
  }

  // Allocate context for decoding the codec
  AVCodecContext * codec_context = avcodec_alloc_context3(codec);
  if (codec_context == NULL) {
    avformat_close_input(&format_context);
    throw std::runtime_error(
        "Could not allocate a decoding context for file: " + filename);
  }

  // Fill the codecContext with parameters of the codec
  if (avcodec_parameters_to_context(
        codec_context, 
        format_context -> streams[audio_stream_index] -> codecpar
        ) != 0) {
    avcodec_close(codec_context);
    avcodec_free_context(&codec_context);
    avformat_close_input(&format_context);
    throw std::runtime_error(
        "Could not set codec context parameters for file: " + filename);
  }

  // Ask for non planar data
  codec_context -> request_sample_fmt =
    av_get_alt_sample_fmt(codec_context -> sample_fmt, 0);

  // Initialize the decoder
  if (avcodec_open2(codec_context, codec, NULL) != 0) {
    avcodec_close(codec_context);
    avcodec_free_context(&codec_context);
    avformat_close_input(&format_context);
    throw std::runtime_error(
        "Could not initialize the decoder for file: " + filename);
  }

  // Allocate a frame
  AVFrame * frame = NULL;
  if ((frame = av_frame_alloc()) == NULL) {
    avcodec_close(codec_context);
    avcodec_free_context(&codec_context);
    avformat_close_input(&format_context);
    throw std::runtime_error(
        "Could not allocate audio frame for file: " + filename);
  }

  // prepare a packet
  AVPacket packet;
  av_init_packet(&packet);

  // END OF FFMPEG SET UP

  // fetch the sample rate
  *sample_rate = codec_context -> sample_rate;

  // Get start and end values in samples
  if (end_seconds == 0)
    end_seconds = (format_context -> duration)/(double)AV_TIME_BASE;
  double start_sample = Units::seconds_to_samples(start_seconds, *sample_rate);
  double end_sample = Units::seconds_to_samples(end_seconds, *sample_rate);

  // Allocate the output vector
  std::vector<std::vector<double>> output(
      codec_context -> channels, 
      std::vector<double>(end_sample - start_sample));

  int read_error;
  size_t sample_offset = 0;

  // Read the file until either nothing is left
  // or we reach desired end of sample
  while (sample_offset < end_sample and
      (read_error = av_read_frame(format_context, &packet)) != AVERROR_EOF) {

    // Are we reading correctly?
    if (read_error != 0) {
      throw std::runtime_error(
          "Read error: " + std::to_string(read_error));
      break;
    }

    // Is this the correct stream?
    if (packet.stream_index != audio_stream_index) {
      // Free the frame buffer and reset
      av_packet_unref(&packet);
      continue;
    }

    // Send the packet to the decoder!
    if (avcodec_send_packet(codec_context, &packet) == 0) {
      // Success!
      // Destroy it :)
      // free the buffers and reset fields
      av_packet_unref(&packet);
    } else  {
      // Failure!
      throw std::runtime_error(
          "Send packet error!");
      break;
    }

    // receive the frame
    while ((read_error = avcodec_receive_frame(codec_context, frame)) == 0) {
      // read the frame
      read_frame(
          codec_context, 
          frame,
          output,
          start_sample,
          end_sample,
          &sample_offset);

      // free buffers and set fields to defaults
      av_frame_unref(frame);
    }

    if (read_error != AVERROR(EAGAIN)) {
      throw std::runtime_error(
          "Receive packet error!");
      break;
    }
  }

  // Drain the decoder

  // Clean up
  // Free all frame data
  av_frame_free(&frame);
  // Close the codec context
  avcodec_close(codec_context);
  // Free the context
  avcodec_free_context(&codec_context);
  avformat_close_input(&format_context);

  return output;
}

void Audio::read_frame(
    const AVCodecContext * codec_context, 
    const AVFrame * frame,
    std::vector<std::vector<double>> & output,
    double start_sample,
    double end_sample,
    size_t * sample_offset) {

  size_t num_samples = frame -> nb_samples;
  size_t num_channels = frame -> channels;

  // for every channel and sample in the frame
  for (size_t channel = 0; channel < num_channels; channel++) {
    for (size_t sample = 0; sample < num_samples; sample++) {

      // If the sample is within our sample range
      if (sample + *sample_offset >= start_sample and
          sample + *sample_offset < end_sample) {

        // get the sample value
        double value = read_sample(codec_context, frame, sample, channel);

        // put it in the output
        output[channel][sample + *sample_offset - start_sample] = value;
      }
    }
  }

  *sample_offset += num_samples;
}

double Audio::read_sample(
    const AVCodecContext * codec_context,
    const AVFrame * frame,
    size_t sample,
    size_t channel) {

  // where it is stored
  uint8_t * frameBuffer;
  // its index
  int sampleLocation;

  // Is it planar?
  // I.E. Is each channel in its own vector
  if (av_sample_fmt_is_planar(codec_context -> sample_fmt)) {
    // yes, we search through the channel
    frameBuffer = frame -> extended_data[channel];
    sampleLocation = sample;
  } else {
    // no, we un-interleave them
    frameBuffer = frame -> extended_data[0];
    sampleLocation = sample * (codec_context -> channels) + channel;
  }

  // the number of bytes in an audio sample
  int sampleBytes = av_get_bytes_per_sample(codec_context -> sample_fmt);

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
      throw std::runtime_error(
          "Sample format is invalid!");
      return 0;
  }

  // if the raw values are padded for uneven bit depths (i.e. 24)
  if (codec_context -> bits_per_raw_sample > 0) {
    rawValue = (rawValue >> (sampleBytes * 8 - codec_context -> bits_per_raw_sample));
  }

  double value;

  switch (codec_context -> sample_fmt) {
    case AV_SAMPLE_FMT_U8:
    case AV_SAMPLE_FMT_U8P:
    case AV_SAMPLE_FMT_S16:
    case AV_SAMPLE_FMT_S16P:
    case AV_SAMPLE_FMT_S32:
    case AV_SAMPLE_FMT_S32P:
      // value / (2 ^ (number of Unsigned Bits) - 1)
      value = rawValue / (double)((1 << (sampleBytes * 8 - 1)) - 1);
      break;
    case AV_SAMPLE_FMT_FLT:
    case AV_SAMPLE_FMT_FLTP:
      value = static_cast<double>(*reinterpret_cast<float *>(&rawValue));
      break;
    case AV_SAMPLE_FMT_DBL:
    case AV_SAMPLE_FMT_DBLP:
      value = *reinterpret_cast<double *>(&rawValue);
      break;
    default:
      throw std::runtime_error(
          "Sample format is invalid!");
      return 0;
  }

  return value;

}
