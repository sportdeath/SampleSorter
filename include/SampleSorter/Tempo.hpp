#ifndef TEMPO_H
#define TEMPO_H

class Tempo {
  public:

  static double tempoToSeconds(double tempo);
  static double tempoToSamples(double tempo, long sampleRate);
  static double secondsToSamples(double seconds, long sampleRate);
  static double samplesToSeconds(double samples, long sampleRate);
  static double tempoToBins(double tempo, long hopSize, long sampleRate);
  static double binsToSamples(long bin, long hopSize);
  static double binsToSeconds(long bin, long hopSize, long sampleRate);

  static double correlationTempo(std::vector<double> offsets,
                                 long hopSize,
                                 long sampleRate);

  static double tempoValue(double tempo,
                           double theOne,
                           std::vector<double> offsets,
                           long hopSize,
                           long sampleRate);


  static std::pair<double, double> fineTuneTempo(double tempo,
                              std::vector<double> offsets,
                              double BPSrange,
                              double BPSstepSize,
                              long hopSize,
                              long sampleRate);
};

#endif
