#ifndef SAMPLE_FILE_H
#define SAMPLE_FILE_H

#include <string>

#include "SampleSorter/AudioSample.hpp"

class SampleFile {
  private:
    AudioSample sample;

    virtual std::vector< std::vector<double> > getWaves() = 0;
    virtual long getSampleRate() const = 0;
  protected:
    std::string filePath;
  public:
    SampleFile(std::string filePath_);
    void process();

    std::string getFilePath();
    std::string getFileName();

    virtual double getSampleLength() const = 0;
    AudioSample * getAudioSample();
};

//extern "C" {
  //SampleFile * SampleFile(char * fileName) {
    //return new SampleFile(fileName);
  //}
  //double getTempo(SampleFile * s) {
    //return s -> getAudioSample -> getTempo();
  //}
//}

#endif
