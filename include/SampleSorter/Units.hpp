#ifndef UNITS_H
#define UNITS_H

#include <vector>
#include <iostream>

class Units {
  public:

    static double tempoToSeconds(double tempo);

    static double tempoToSamples(double tempo, long sampleRate);

    static double secondsToSamples(double seconds, long sampleRate);

    static double secondsToTempo(double seconds);

    static double samplesToSeconds(double samples, long sampleRate);

    static double tempoToBins(double tempo, long hopSize, long sampleRate);

    static double binsToSamples(double bin, long hopSize);

    static double binsToSeconds(double bin, long hopSize, long sampleRate);

    static double binsToTempo(double bin, long hopSize, long sampleRate);

    static double samplesToBins(double samples, long hopSize);

    static double secondsToBins(double seconds, long hopSize, long sampleRate);
};

#endif
