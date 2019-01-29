#pragma once

namespace sample_sorter {

class Units {
  public:

    static double tempo_to_seconds(double tempo) {
      return 1/tempo;
    }

    static double seconds_to_tempo(double seconds) {
      return 1/seconds;
    }

    static double seconds_to_samples(double seconds, double sample_rate) {
      return seconds * sample_rate;
    }

    static double samples_to_seconds(double samples, double sample_rate) {
      return samples/sample_rate;
    }

    static double tempo_to_samples(double tempo, double sample_rate) {
      return seconds_to_samples(tempo_to_seconds(tempo), sample_rate);
    }

    static double samples_to_tempo(double samples, double sample_rate) {
      return seconds_to_tempo(samples_to_seconds(samples, sample_rate));
    }
};

}
