#ifndef SAMPLE_FILE_H
#define SAMPLE_FILE_H

#include <string>

#include "SampleSorter/AudioSample.hpp"

class SampleFile {
  private:
    virtual std::vector< std::vector<double> > extractAudio(long * sampleRate) = 0;
    virtual bool readMetaData() = 0;
  protected:
    std::string filePath;
    AudioSample sample;
  public:
    SampleFile(std::string filePath_);
    bool process(bool forceReprocess);

    std::string getFilePath();
    std::string getFileName();

    virtual double getSampleSeconds() const = 0;
    AudioSample * getAudioSample();
};

#endif
