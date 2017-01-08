#ifndef TEMPO_FUNCTION_H
#define TEMPO_FUNCTION_H

#include <vector>
#include <cmath>
#include <random>

class TempoFunction {
  private:
    double avgTempo;
    double totalSeconds;
    std::vector<double> aCoefficients;
    std::vector<double> bCoefficients;
  public:
    TempoFunction();
    double getTotalSeconds() const;
    double getTheOne() const;

    TempoFunction(double avgTempo_,
                  double theOne,
                  double totalSeconds_);

    TempoFunction(
        std::vector<std::vector<double> > & audio,
        long sampleRate,
        size_t degrees, 
        double percentageError,
        double percentageStep);

    void fineTuneAvgTempo(
        double guessTempo,
        double percentageError,
        double percentageStep,
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate);

    void gradientDescent(
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate
        );

    void fineTuneTheOne(
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate
        );

    double getAvgTempo() const;

    size_t getDegree() const;

    double getBeatNum(double time) const;

    // The derivative
    double getTempo(double time) const;

    double distanceFromBeat(double time) const;

    double getValue(
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate
        ) const;

    void getBinGradient(
        size_t bin,
        std::vector<double> & aNabla,
        std::vector<double> & bNabla,
        const std::vector<double> & onsets,
        long hopSize,
        long sampleRate
        ) const;

    double getKthBeat(double k, double accuracy = 0.005) const;

    static double correlationTempo(
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
