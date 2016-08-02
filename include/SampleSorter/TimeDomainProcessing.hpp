#ifndef TIME_DOMAIN_PROCESSING_H
#define TIME_DOMAIN_PROCESSING_H

class TimeDomainProcessing {
  public:
    static double getEnergy(const std::vector<double> & signal);
    static double getAverageEnergyPerSample(const std::vector<double> & signal);
    static double getAverageEnergyPerSecond(const std::vector<double> & signal,
                                            long sampleRate);
    static void normalizeEnergy(std::vector< std::vector<double> > & output, 
                                const std::vector< std::vector<double> > & audio,
                                long sampleRate);
};

#endif
