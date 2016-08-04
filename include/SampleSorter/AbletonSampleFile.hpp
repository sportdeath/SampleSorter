#ifndef ABLETON_SAMPLE_FILE_H
#define ABLETON_SAMPLE_FILE_H

#include <vector>

#include <tinyxml2.h>

#include <sndfile.hh>

#include "SampleSorter/SampleFile.hpp"

class AbletonSampleFile : public SampleFile {
  private:
    tinyxml2::XMLDocument doc;

    std::string referenceFilePath;
    double startSeconds;
    double endSeconds;

    SndfileHandle audioFile;

    void getDoc();
    void readDoc();

    std::vector< std::vector<double> > getWaves();
    long getSampleRate() const;

  public:
    AbletonSampleFile(std::string filePath);

    AbletonSampleFile(const AbletonSampleFile & other);

    std::string getReferenceFilePath() const;
    std::string getReferenceFileName() const;

    double getSampleLength() const;
};

#endif
