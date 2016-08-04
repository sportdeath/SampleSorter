#ifndef AUDIO_SAMPLE_H
#define AUDIO_SAMPLE_H

#include <vector>

#include "SampleSorter/Octave.hpp"

class AudioSample {
  private:
    long tuningCents;
    double tempo;
    double theOne;
    std::vector<Octave> chords;

    void tune(std::vector<std::vector<double> > audio, long sampleRate);
    void findBeat(std::vector<std::vector<double> > audio, long sampleRate);
    void findChords(std::vector<std::vector<double> > audio, long sampleRate);
  public:
    AudioSample();
    AudioSample(std::vector<std::vector<double> > audio, long sampleRate);

    long getTuningCents() const;
    double getTuningCentsFreqRatio() const;
    double getBeatRaw() const;
    double getBeatWithTuning() const;
    double getTheOneRaw() const;
    double getTheOneWithTuning() const;
    std::vector<Octave> getChords() const;
};

#endif
