#ifndef AUDIO_SAMPLE_H
#define AUDIO_SAMPLE_H

#include <vector>

#include "SampleSorter/Octave.hpp"
#include "SampleSorter/TempoFunction.hpp"

class AudioSample {
  private:
    long tuningCents;
    TempoFunction tempo;
    std::vector<Octave> chords;

    void tune(std::vector<std::vector<double> > & audio, long sampleRate);
    void findBeat(std::vector<std::vector<double> > & audio, long sampleRate);
    void findChords(std::vector<std::vector<double> > & audio, long sampleRate);
  public:
    AudioSample();
    AudioSample(long tuningCents_,
                double rawBeat,
                double theOne,
                double totalSeconds,
                std::vector<Octave> chords);
    AudioSample(std::vector<std::vector<double> > & audio, long sampleRate);

    double getTotalSeconds() const;
    long getTuningCents() const;
    double getTuningCentsFreqRatio() const;
    double getBeatRaw() const;
    double getBeatWithTuning() const;
    double getTheOneRaw() const;
    double getTheOneWithTuning() const;
    std::vector<Octave> getChords() const;
};

#endif
