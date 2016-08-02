#ifndef ABLETON_SAMPLE_H
#define ABLETON_SAMPLE_H

#include <vector>

#include <tinyxml2.h>

#include "SampleSorter/Sample.hpp"

class AbletonSample : public Sample {
  private:
    std::string name;

    bool docExists = false;
    tinyxml2::XMLDocument doc;

    bool wavesExist = false;
    std::vector< std::vector<double> > waves;
    long sampleRate;

    tinyxml2::XMLDocument * getDoc();

    std::vector< std::vector<double> > getWaves();
    long getSampleRate();

  public:

    std::string getName() const;

    AbletonSample(std::string file);

    AbletonSample(const AbletonSample & other);
};

#endif
