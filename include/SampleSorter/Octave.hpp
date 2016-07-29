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

    static const std::vector<double> harmonicSeries;

    static const std::vector<std::string> noteNames;

    static const Octave unison;
    static const Octave minorSecond;
    static const Octave majorSecond;
    static const Octave minorThird;
    static const Octave majorThird;
    static const Octave perfectFifth;
    
    static const Octave majorTriad;
    static const Octave augmentedTriad;
    static const Octave minorTriad;
    static const Octave diminishedTriad;

    static const std::vector<Octave> chords;

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

    void addPeak(double peakFreq, double peakValue, double baseOffset);

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

    void add(const Octave & other);
    void normalize();

    // 
    double difference (const Octave & other);
    double similarity (const Octave & other);
    double dissonance ();

    static Octave mean(
        std::vector<Octave>::iterator begin,
        std::vector<Octave>::iterator end);
    static double selfSimilarity(
        std::vector<Octave>::iterator begin,
        std::vector<Octave>::iterator end);
    static double dissonance(
        std::vector<Octave>::iterator begin,
        std::vector<Octave>::iterator end);

    std::string Octave::mostSimilarChord();
};

#endif
