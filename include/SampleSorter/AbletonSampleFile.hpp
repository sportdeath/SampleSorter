#ifndef ABLETON_SAMPLE_FILE_H
#define ABLETON_SAMPLE_FILE_H

#include <vector>

#include <tinyxml2.h>

#include "SampleSorter/SampleFile.hpp"

//class AbletonNode {
  //private:
  //public:
    //// constructor with list
    //AbletonNode
    

//}

class AbletonSampleFile : public SampleFile {
  private:
    tinyxml2::XMLDocument doc;

    std::string userLibrary;

    std::string referenceFilePath;
    double startSeconds;
    double endSeconds;

    // opens file
    void getDoc();

    // parses file
    // returns true iff preprocessed
    bool readDoc();

    // writes file
    void setDoc();

    virtual std::vector< std::vector<double> > extractAudio(long * sampleRate);

    // opens and parses file
    virtual bool readMetaData();

    // Constants
    //static const std::string VALUE = "Value";

    // Document reading functions
    tinyxml2::XMLElement * getNode(std::initializer_list<std::string> nodes);
    tinyxml2::XMLElement * getAudioNode();
    tinyxml2::XMLElement * getLoopNode();
    tinyxml2::XMLElement * getSortDataNode();

    // Direct ones
    tinyxml2::XMLElement * getLoopStartNode();
    tinyxml2::XMLElement * getLoopEndNode();
    tinyxml2::XMLElement * getHiddenLoopStartNode();
    tinyxml2::XMLElement * getHiddenLoopEndNode();
    tinyxml2::XMLElement * getNameNode();
    tinyxml2::XMLElement * getPitchFineNode();
    
    // various fetchers
    std::string fetchReferenceFilePath();

  public:
    AbletonSampleFile(std::string filePath, std::string userLibrary);

    AbletonSampleFile(const AbletonSampleFile & other);

    std::string getReferenceFilePath() const;
    std::string getReferenceFileName() const;

    virtual double getSampleSeconds() const;

    void writeToFile();
};

#endif
