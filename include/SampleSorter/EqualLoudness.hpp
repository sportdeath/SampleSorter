#ifndef EQUAL_LOUDNESS_H
#define EQUAL_LOUDNESS_H

class EqualLoudness {
  public:
    static double AWeightingAmp(double freq);
    static double CWeightingAmp(double freq);
    static void filter(
        std::vector<std::vector<double> > & output, 
        std::vector<std::vector<double> > & audio, 
        long sampleRate);

};

#endif
