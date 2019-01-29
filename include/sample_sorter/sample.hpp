#pragma once

#include <cmath>
#include <vector>

#include "sample_sorter/octave.hpp"
#include "sample_sorter/beat.hpp"

namespace sample_sorter {

/**
 * A high level representation of an audio sample
 * that is removed from the actual audio file.
 *
 * This contains information about the audio sample's
 * tempo, tuning and harmony which can be used to classify
 * the samples.
 */
class Sample {

  private:

    // The length in seconds of the sample
    double total_seconds_;

    // The sample's tuning relative to A440
    double tuning_cents_;

    // A representation of the sample's tempo and phase
    Beat beat;

    // A vector of octaves for each beat of the sample
    std::vector<Octave> chords_;

    // Tune the sample
    void tune(std::vector<std::vector<double> > & audio, double sample_rate);

    // Find the tempo of the sample
    void find_tempo(std::vector<std::vector<double> > & audio, double sample_rate);

    // Find a harmonic representation of the sample
    void find_chords(std::vector<std::vector<double> > & audio, double sample_rate);

  public:
    // Empty constructor
    Sample() {};

    // Fill in the parameters with audio
    Sample(std::vector<std::vector<double> > & audio, double sample_rate);

    // Used for loading preexisting data
    Sample(
        double _total_seconds,
        double _tuning_cents,
        double tempo,
        double the_one,
        std::vector<Octave> _chords)
      : total_seconds_(_total_seconds),
        tuning_cents_(_tuning_cents),
        beat(tempo, the_one),
        chords_(_chords) {}

    // Publicize class members
    double total_seconds() const {return total_seconds_;}
    double tuning_cents() const {return tuning_cents_;}
    const std::vector<Octave> & chords() const {return chords_;}
    double tuning_cents_freq_ratio() const {
      return std::pow(2., tuning_cents()/1200.);}
    double tempo() const {return beat.tempo();}
    double tempo_with_tuning() const {
      return tempo() * tuning_cents_freq_ratio();}
    double the_one() const {return beat.the_one();}
    double the_one_with_tuning() const {
      return the_one()/tuning_cents_freq_ratio();}
};

}
