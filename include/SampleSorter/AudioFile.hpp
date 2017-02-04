#ifndef AUDIO_FILE_H
#define AUDIO_FILE_H

#include <vector>
#include <string>

#include <libavcodec/avcodec.h>

class AudioFile {
  public:
    static std::vector<std::vector<double> > read(
        std::string fileName,
        double startSeconds,
        double endSeconds,
        long * sampleRate);

    static void readFrame(
        const AVCodecContext * codecContext,
        const AVFrame * frame,
        std::vector<std::vector<double> > & output,
        long startSample,
        long endSample,
        long * sampleOffset);

    static double getSample(
        const AVCodecContext * codecContext,
        const AVFrame * frame,
        int sample,
        int channel);
};

#endif
