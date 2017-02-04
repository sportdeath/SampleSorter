#ifndef AUDIO_SAMPLE_H
#define AUDIO_SAMPLE_H

#include <vector>

#include "SampleSorter/Octave.hpp"
#include "SampleSorter/Tempo.hpp"

class AudioSample {
  private:
    long tuningCents;
    Tempo tempo;
    std::vector<Octave> chords;
    double totalSeconds;
    long sampleRate;

    void tune(std::vector<std::vector<double> > & audio);
    void findBeat(std::vector<std::vector<double> > & audio);
    void findChords(std::vector<std::vector<double> > & audio);

    static const double TEMPO_PERCENTAGE_ERROR;
    static const double TEMPO_STEPS;
    static const double TEMPO_ONE_STEPS;
  public:
    AudioSample();
    AudioSample(long tuningCents_,
                double rawBeat,
                double theOne,
                double totalSeconds,
                long sampleRate,
                std::vector<Octave> chords);
    AudioSample(std::vector<std::vector<double> > & audio, long _sampleRate);

    double getTotalSeconds() const;
    long getTuningCents() const;
    double getTuningCentsFreqRatio() const;
    double getBeatRaw() const;
    double getBeatWithTuning() const;
    double getTheOneRaw() const;
    double getTheOneWithTuning() const;
    std::vector<Octave> getChords() const;
    long getSampleRate() const;
    double getLastBeatSeconds() const;
};

#endif
