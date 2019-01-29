#pragma once

#include <vector>
#include <cmath>
#include <random>

namespace sample_sorter {

class Beat {
  private:
    /**
     * The average tempo of the
     * measured in beats per second
     */
    double tempo_;

    /**
     * The location of the starting beat
     * in the sample measured in seconds
     */
    double the_one_;

  public:
    Beat() {}
    Beat(
        double tempo, 
        double the_one)
      : tempo_(tempo),
        the_one_(the_one) {}

    /**
     * Takes an audio file and measures
     * the tempo and the one.
     */
    Beat(
        const std::vector<std::vector<double>> & audio,
        double sample_rate);

    std::vector<double> onset_energy(
        const std::vector<std::vector<double>> & audio,
        double sample_rate) const;

    double the_one() const {return the_one_;}
    double tempo() const {return tempo_;}
};

}
