#include <iostream>

#include "SampleSorter/SampleFile.hpp"
#include "SampleSorter/AbletonSampleFile.hpp"
#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/Octave.hpp"
#include "Plotting/Plotting.hpp"

extern "C" {
  SampleFile * NewAbletonSampleFile(const char * filePath, const char * userLibrary) {
    AbletonSampleFile * s = new AbletonSampleFile(filePath, userLibrary);
    return s;
  }
  bool process(SampleFile * s, bool forceReprocess) {
    return s -> process(forceReprocess);
  }
  const char * getFileName(SampleFile * s) {
    return s -> getFileName().c_str();
  }
  long getTuningCents(SampleFile * s) {
    return s -> getAudioSample() -> getTuningCents();
  }
  short getFundemental(SampleFile * s) {
    return s -> getAudioSample() -> getFundemental();
  }
  long getTheOneWithTuning(SampleFile * s) {
    return s -> getAudioSample() -> getTheOneWithTuning();
  }
  double getBeatWithTuning(SampleFile * s) {
    return s -> getAudioSample() -> getBeatWithTuning();
  }
  size_t getNumChords(SampleFile * s) {
    return s -> getAudioSample() -> getChords().size();
  }
  double ** getChords(SampleFile * s) {
    std::vector<Octave> chords = s -> getAudioSample() -> getChords();
    double ** chordsArray = new double *[chords.size()];
    for (long i = 0; i < chords.size(); i++) {
      chordsArray[i] = new double[12];
      for (long j = 0; j < 12; j++) {
        chordsArray[i][j] = chords[i].getSpectrogram()[j];
      }
    }
    return chordsArray;
  }
  void writeToFile(AbletonSampleFile * s) {
    s -> writeToFile();
  }
  void deleteChords(double ** chords, size_t chordsSize) {
    for (size_t i = 0; i < chordsSize; i++) {
      delete chords[i];
    }
    delete chords;
  }
  void deleteAbletonSampleFile(AbletonSampleFile * s) {
    delete s;
  }
}
