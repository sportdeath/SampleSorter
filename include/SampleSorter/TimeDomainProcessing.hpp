#ifndef TIME_DOMAIN_PROCESSING_H
#define TIME_DOMAIN_PROCESSING_H

class TimeDomainProcessing {
  public:
    static double getEnergy(const std::vector<double> & signal);
    static double getAverageEnergyPerSample(const std::vector<double> & signal);
    static double getAverageEnergyPerSecond(const std::vector<double> & signal,
                                            long sampleRate);
    static double getAverageEnergyPerBeat(const std::vector<double> & signal,
                                          double tempo,
                                          long sampleRate);
    static void normalizeByEnergy(std::vector< std::vector<double> > & output, 
                                const std::vector< std::vector<double> > & audio,
                                double energy);
    static void unitEnergyPerSecond(std::vector< std::vector<double> > & output, 
                                const std::vector< std::vector<double> > & audio,
                                long sampleRate);
    static void unitEnergyPerBeat(std::vector< std::vector<double> > & output, 
                                const std::vector< std::vector<double> > & audio,
                                double tempo,
                                long sampleRate);
};

#endif
