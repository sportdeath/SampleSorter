#pragma once

#include <vector>
#include <string>

#include <libavcodec/avcodec.h>

namespace sample_sorter {

class Audio {
  private:
    static void read_frame(
        const AVCodecContext * codec_context,
        const AVFrame * frame,
        std::vector<std::vector<double>> & output,
        double start_sample,
        double end_sample,
        size_t * sampleOffset);

    static double read_sample(
        const AVCodecContext * codec_context,
        const AVFrame * frame,
        size_t sample,
        size_t channel);

  public:

    static std::vector<std::vector<double>> read(
        std::string filename,
        double * sample_rate,
        double start_seconds=0,
        double end_seconds=0);

    static void write(
        const std::vector<std::vector<double>> & audio,
        std::string filename,
        double sample_rate);
};

}
