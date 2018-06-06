#include <vector>
#include <iostream>
#include <ciso646>

#include <boost/filesystem.hpp>

#include "SampleSorter/SampleFile.hpp"
#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/ProcessingException.hpp"

SampleFile::SampleFile(std::string filePath_) {
  filePath = filePath_;
}

bool SampleFile::process() {
  try {
    if (not readMetaData()) {
      long sampleRate;
      std::vector< std::vector<double> > waves = extractAudio(&sampleRate);
      sample = AudioSample(waves, sampleRate);
    }
    return true;
  } catch (ProcessingException & e) {
    std::cout << e.getMessage() << std::endl;
    return false;
  }
}

std::string SampleFile::getFilePath() {
  return filePath;
}

std::string SampleFile::getFileName() {
  return boost::filesystem::path(filePath).stem().make_preferred().string();
}

AudioSample * SampleFile::getAudioSample() {
  return &sample;
}
