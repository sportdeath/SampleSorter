#include <vector>

#include <tinyxml2.h>

#include "SampleSorter/Sample.hpp"

class AbletonSample : public Sample {
  private:
    bool docExists;
    tinyxml2::XMLDocument doc;
    bool wavesExist;
    std::vector< std::vector<double> > waves;
    long sampleRate;

    tinyxml2::XMLDocument * getDoc();

    std::vector< std::vector<double> > getWaves();
    long getSampleRate();

  public:
    AbletonSample(std::string file);
};
