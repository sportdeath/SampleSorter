#ifndef ABLETON_SAMPLE_FILE_H
#define ABLETON_SAMPLE_FILE_H

#include <vector>

#include <tinyxml2.h>

#include "SampleSorter/SampleFile.hpp"

class AbletonSampleFile : public SampleFile {
  private:
    tinyxml2::XMLDocument doc;

    std::string referenceFilePath;
    double startSeconds;
    double endSeconds;

    void getDoc();
    // returns true iff preprocessed
    bool readDoc();

    std::vector< std::vector<double> > extractAudio(long * sampleRate);

    bool readMetaData();

  public:
    AbletonSampleFile(std::string filePath);

    AbletonSampleFile(const AbletonSampleFile & other);

    tinyxml2::XMLElement * getAudioNode();
    tinyxml2::XMLElement * getLoopNode();

    std::string getReferenceFilePath() const;
    std::string getReferenceFileName() const;

    double getSampleLength() const;

    void writeToFile();
};

#endif
