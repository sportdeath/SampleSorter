#ifdef _WIN32
#  define EXPORT __declspec(dllexport)
#else
#  define EXPORT
#endif

#include <iostream>

#include "SampleSorter/SampleFile.hpp"
#include "SampleSorter/AbletonSampleFile.hpp"
#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/Octave.hpp"
#include "Plotting/Plotting.hpp"

extern "C" {
  EXPORT AbletonSampleFile * NewAbletonSampleFile(const char * filePath, const char * userLibrary, bool forceReprocess) {
    AbletonSampleFile * s = new AbletonSampleFile(filePath, userLibrary, forceReprocess);
    return s;
  }
  EXPORT bool process(SampleFile * s) {
    return s -> process();
  }
  EXPORT const char * getFileName(SampleFile * s) {
    return s -> getFileName().c_str();
  }
  EXPORT const char * getAudioPath(AbletonSampleFile * s) {
    return s -> getReferenceFilePath().c_str();
  }
  EXPORT long getTuningCents(SampleFile * s) {
    return s -> getAudioSample() -> getTuningCents();
  }
  EXPORT short getFundemental(SampleFile * s) {
    return s -> getAudioSample() -> getFundemental();
  }
  EXPORT long getTheOneWithTuning(SampleFile * s) {
    return s -> getAudioSample() -> getTheOneWithTuning();
  }
  EXPORT double getBeatWithTuning(SampleFile * s) {
    return s -> getAudioSample() -> getBeatWithTuning();
  }
  EXPORT size_t getNumChords(SampleFile * s) {
    return s -> getAudioSample() -> getChords().size();
  }
  EXPORT double ** getChords(SampleFile * s) {
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
  EXPORT double * getOctave(SampleFile * s) {
    Octave octave = s -> getAudioSample() -> getOctave();
    double * octaveArray = new double [12];
    for (long i = 0; i < 12; i++) {
      octaveArray[i] = octave.getSpectrogram()[i];
    }
    return octaveArray;
  }
  EXPORT void writeToFile(AbletonSampleFile * s) {
    s -> writeToFile();
  }
  EXPORT void deleteChords(double ** chords, size_t chordsSize) {
    for (size_t i = 0; i < chordsSize; i++) {
      delete chords[i];
    }
    delete chords;
  }
  EXPORT void deleteOctave(double * octave) {
    delete octave;
  }
  EXPORT void deleteAbletonSampleFile(AbletonSampleFile * s) {
    delete s;
  }
}
