#ifndef OCTAVE_H
#define OCTAVE_H

#include <cmath>
#include <vector>

#include <fftw3.h>

class Octave {
  private:
    // A vector with magnitude 1
    std::vector<double> spectrogram;
  public:
    static const int BASE_FREQ = 440;

    static const int SEMITONES_PER_OCTAVE = 12;
    static const int CENTS_PER_SEMITONE = 100;
    static const int CENTS_PER_OCTAVE = SEMITONES_PER_OCTAVE * CENTS_PER_SEMITONE;
    double CENTS_PER_BIN;
    double BINS_PER_SEMITONE;

    double magnitude = 1;

    Octave(std::vector< std::vector<double> > audio, 
           long numBins,
           long sampleRate,
           long tuningCents = 0);

    Octave(fftw_complex * fft, 
           long fftSize,
           long numBins,
           long sampleRate,
           long tuningCents = 0);

    Octave(long numBins = 12);
    Octave(std::vector<double> spec);

    long getBins() const;

    double getBinsPerSemitone() const;
    double getCentsPerBin() const;

    void addPeak(double peakFreq, double peakValue, double baseOffset);

    void rotate(long cents);

    /**
     * Returns the proportion of spectral energy
     * outside of bins.
     */
    double tuningValue() const;

    /**
     * Tunes the octave by rotating it +-50 cents
     * and choosing the rotation that minimizes the
     * tuning value. The rotation is returned.
     */
    long tune();

    void plot(std::string title, std::string xaxis, std::string yaxis, bool histogram=false) const;

    void add(Octave & output, const Octave & that) const;

    std::vector<double> getSpectrogram() const;

    long getMax() const;
};

#endif
