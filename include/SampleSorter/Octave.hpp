#include <cmath>
#include <vector>

class Octave {
  private:
    std::vector<double> spectrogram;
  public:
    static const int BASE_FREQ = 440;
    double BASE_OFFSET;

    static const int SEMITONES_PER_OCTAVE = 12;
    static const int CENTS_PER_SEMITONE = 100;
    static const int CENTS_PER_OCTAVE = SEMITONES_PER_OCTAVE * CENTS_PER_SEMITONE;
    double CENTS_PER_BIN;
    double BINS_PER_SEMITONE;

    Octave(std::vector< std::vector<double> > audio, long sampleRate, long numBins);

    void addPeak(double peakFreq, double peakValue);

    void rotate(long cents);

    /**
     * Returns the proportion of spectral energy
     * outside of bins.
     */
    double tuningValue();

    /**
     * Tunes the octave by rotating it +-50 cents
     * and choosing the rotation that minimizes the
     * tuning value. The rotation is returned.
     */
    long tune();

    void plot();
};
