#pragma once

#include <vector>

#include "SampleSorter/Octave.hpp"
#include "SampleSorter/Tempo.hpp"

/**
 * A high level representation of an audio sample
 * that is removed from the actual audio file.
 *
 * This contains information about the audio sample's
 * tempo, tuning and harmony which can be used to classify
 * the samples.
 */
class AudioSample {

  private:

    // The length in seconds of the sample
    double totalSeconds;

    // The sample's tuning relative to A440
    long tuningCents;

    // The sample's tempo
    Tempo tempo;

    // A vector of octaves for each beat of the sample
    std::vector<Octave> chords;

    // Tune the sample
    void tune(std::vector<std::vector<double> > & audio, long sampleRate);

    // Find the tempo of the sample
    void findTempo(std::vector<std::vector<double> > & audio, long sampleRate);

    // Find a harmonic representation of the sample
    void findChords(std::vector<std::vector<double> > & audio, long sampleRate);

  public:
    // Empty constructor
    AudioSample() {};

    // Fill in the parameters with audio
    AudioSample(std::vector<std::vector<double> > & audio, long sampleRate);

    // Used for loading preexisting data
    AudioSample(
        long tuningCents_,
        double tempo_,
        double theOne,
        double totalSeconds_,
        std::vector<Octave> chords_)
      : tuningCents(tuningCents_),
        tempo(tempo_, theOne),
        totalSeconds(totalSeconds_),
        chords(chords_) {}

    // Getters
    double getTotalSeconds() const {return totalSeconds;}
    long getTuningCents() const {return tuningCents;}
    double getTuningCentsFreqRatio() const {
      return std::pow(2., tuningCents/1200.);}
    double getTempo() const {return tempo.getTempo();}
    double getTempoWithTuning() const {
      return tempo.getTempo() * getTuningCentsFreqRatio();}
    double getTheOne() const {return tempo.getTheOne();}
    double getTheOneWithTuning() const {
      return tempo.getTheOne()/getTuningCentsFreqRatio();}
    const std::vector<Octave> & getChords() const {return chords;}
};
