#include <vector>
#include <iostream>

#include <boost/filesystem.hpp>

#include "SampleSorter/SampleFile.hpp"
#include "SampleSorter/AudioSample.hpp"

SampleFile::SampleFile(std::string filePath_) {
  filePath = filePath_;
}

void SampleFile::process() {
  sample = AudioSample(getWaves(), getSampleRate());
}

std::string SampleFile::getFilePath() {
  return filePath;
}

std::string SampleFile::getFileName() {
  return boost::filesystem::path(filePath).stem().native();
}

AudioSample * SampleFile::getAudioSample() {
  return &sample;
}

//void SampleFile::writeToMIDI(std::string fileName) {
  //// Find max
  //double max = 0;
  //for (long i = 0; i < chords.size(); i++) {
    //for (long bin = 0; bin < 12; bin ++) {
      //getSpectrogram;
      //max = std::max(max, spec[bin]);
    //}
  //}

  //ratio = 127./max;

  //// open midi file
  //// set tempo
  //// add track
  //// for every chord
  //// for every bin
  //// add an event
//}




//bool Sample::isCompatible(const Sample & other) {
  //// get larger tempo
  //double larger, smaller;
  //larger = std::max(getBeatWithTuning(), other.getBeatWithTuning());
  //smaller = std::min(getBeatWithTuning(), other.getBeatWithTuning());
  //long ratio = std::round(larger/smaller);
  //double tuningSteps = 12 * std::log2(larger/(ratio * smaller));
  ////double integerDeviation = std::fmod(tuningSteps, 1.);
  //double integerDeviation = tuningSteps - std::round(tuningSteps);

  //if (std::abs(integerDeviation) < 0.05) {
    //std::cout << "'" << getName() << "' is compatible with '" << other.getName() <<"'" << std::endl;
    //std::cout << "Tempos: " << 60*getBeatWithTuning() <<", " << 60*other.getBeatWithTuning() << std::endl;
    //std::cout << "Tuning cents: " << tuningCents <<", " << other.tuningCents << std::endl;
    //std::cout << "Tuning ratio: " << std::round(tuningSteps) << " steps and " << integerDeviation * 100 << " cents" << std::endl << std::endl;
  //}

  //return integerDeviation < 0.05;
//}

  // find whether that ratio is an even number of cents

//bool Sample::isSimilar(const & Sample other) {
  //// coarse tune so their tempos are equal
  //// or are integer ratios of each other
  //// get similarity of all octaves.
  //// Look for large chains
//}
