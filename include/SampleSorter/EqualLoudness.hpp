#ifndef EQUAL_LOUDNESS_H
#define EQUAL_LOUDNESS_H

class EqualLoudness {
  public:
    static double AWeightingAmp(double freq);
    static double CWeightingAmp(double freq);
    static std::vector<std::vector<double> > filter(
        std::vector<std::vector<double> > audio, 
        long sampleRate);

};

#endif
