#pragma once

#include <cmath>
#include <vector>

namespace sample_sorter {

class Octave {
  private:
    std::vector<double> spectrogram;

  public:
    Octave(
        std::vector< std::vector<double> > audio, 
        size_t num_bins,
        double sample_rate,
        double tuning_cents=0);

    // Default initialization
    Octave() {};
    Octave(std::vector<double> spectrogram_) 
        : spectrogram(spectrogram_) {}

    /**
     * Tunes the octave by rotating it
     * and choosing the rotation that minimizes the
     * tuning value. The rotation is returned.
     */
    double tune();

    size_t bins() const {return spectrogram.size();}

    static const int SEMITONES_PER_OCTAVE = 12;
    static const int CENTS_PER_SEMITONE = 100;
    static const int CENTS_PER_OCTAVE = SEMITONES_PER_OCTAVE * CENTS_PER_SEMITONE;
    static const int BASE_FREQ = 440;
};

}
