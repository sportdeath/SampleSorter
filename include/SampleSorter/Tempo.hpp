#ifndef TEMPO_FUNCTION_H
#define TEMPO_FUNCTION_H

#include <vector>
#include <cmath>
#include <random>

class Tempo {
  private:
    /**
     * The average tempo of the
     * sample
     */
    double tempo;
    /**
     * The starting beat of the
     * sample measured in bins
     */
    double theOneBin;
  public:
    /**
     * The tempo.
     */
    double getTempo() const;
    /**
     * Returns the starting
     * beat in bins
     */
    double getTheOneBin() const;
    /**
     * Returns the starting beat
     * in seconds
     */
    double getTheOne(long hopSize, long sampleRate) const;

    Tempo();
    Tempo(double tempo_, double theOneBin_);

    /**
     * Takes an audio file
     * and measures the
     * tempo and the one
     */
    Tempo(
        const std::vector<std::vector<double> > & audio,
        long sampleRate,
        double percentageError,
        int tempoSteps,
        int oneSteps,
        );

    /**
     * For a given guess tempo,
     * this finds the tempo whose length
     * in seconds is within the percentage
     * error bounds that maximizes the energy
     * around beats. This energy is measured
     * by the getValue function.
     */
    void fineTuneTempo(
        const double percentageError,
        const int steps,
        const std::vector<double> & onsets,
        const long hopSize,
        const long sampleRate);

    /**
     * Fine tunes the one so that most of the energy
     * lies after each beat.
     */
    void fineTuneTheOne(
        const std::vector<double> & onsets,
        const int steps,
        const long hopSize,
        const long sampleRate
        );


    /**
     * Computes the distance from a beat.
     * If bidirectional, the distance is,
     * from either the previous or next beat,
     * whichever is closer.
     * If not bidirectional, the distance is
     * the distance from the previous beat.
     */
    double distanceFromBeat(double bin, bool bidirectional, long hopSize, long sampleRate) const;

    double getValue(
        const std::vector<double> & onsets,
        bool bidirectional,
        long hopSize,
        long sampleRate
        ) const;

    void findCorrelationTempo(
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate);

    void plotBeats(double resolution) const;
    void plotTempo(double resolution) const;
    void plotOnsetsWithBeats(
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate) const;
    
};

#endif
