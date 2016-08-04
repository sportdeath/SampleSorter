#include <vector>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/TimeDomainProcessing.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/EqualLoudness.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/Tempo.hpp"

AudioSample::AudioSample() {
  tuningCents = 0;
  tempo = 1;
  theOne = 0;
}

AudioSample::AudioSample(
    std::vector<std::vector<double> > audio, 
    long sampleRate) {

  tune(audio, sampleRate);
  findBeat(audio, sampleRate);
  findChords(audio, sampleRate);
}

void AudioSample::tune(std::vector<std::vector<double> > audio, long sampleRate) {
  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);

  Octave oct(filteredAudio, 1200., sampleRate);
  tuningCents = oct.tune();
}

void AudioSample::findBeat(std::vector<std::vector<double> > audio, long sampleRate) {
  long hopSize = 1024;
  long windowRatio = 2;

  std::vector<double> onsets = 
    SpectralProcessing::onsetEnergy(audio, hopSize, windowRatio);
  tempo = Tempo::correlationTempo(onsets, hopSize, sampleRate);

  std::pair<double, double> tempoOne = 
    Tempo::fineTuneTempo(tempo, onsets, 10/60., 0.05/60., hopSize, sampleRate);

  tempo = tempoOne.first;
  theOne = tempoOne.second;
}

void AudioSample::findChords(std::vector<std::vector<double> > audio, long sampleRate) {
  long windowSize = Tempo::tempoToSamples(tempo, sampleRate);
  std::vector<double> window(windowSize);

  // set up fourier transform
  long fftSize = windowSize/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
                                           window.data(),
                                           fft,
                                           FFTW_ESTIMATE);

  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);
  TimeDomainProcessing::unitEnergyPerBeat(
      filteredAudio, filteredAudio, tempo, sampleRate);

  long theOneSamples = Tempo::secondsToSamples(theOne, sampleRate);
  long maxWindow = (filteredAudio[0].size() - theOneSamples)/windowSize;

  chords.resize(maxWindow);

  for (long hop = 0; hop < maxWindow; hop++) {
    for (long channel = 0; channel < filteredAudio.size(); channel ++) {
      for (long i = 0; i < windowSize; i++) {
        window[i] = filteredAudio[channel][theOneSamples + i + hop * windowSize];
        window[i] = window[i] * SpectralProcessing::hammingWindow(i, windowSize);
      }

      fftw_execute(fftPlan);
      Octave channelOctave(fft, fftSize, 12, sampleRate, tuningCents);

      chords[hop].add(chords[hop], channelOctave);
    }
  }

  fftw_free(fft);
  fftw_destroy_plan(fftPlan);
}

long AudioSample::getTuningCents() const {
  return tuningCents;
}
double AudioSample::getTuningCentsFreqRatio() const {
  return std::pow(2., tuningCents/1200.);
}
double AudioSample::getBeatRaw() const {
  return tempo;
}
double AudioSample::getBeatWithTuning() const {
  return tempo * getTuningCentsFreqRatio();
}
double AudioSample::getTheOneRaw() const {
  return theOne;
}
double AudioSample::getTheOneWithTuning() const {
  return theOne/getTuningCentsFreqRatio();
}
std::vector<Octave> AudioSample::getChords() const {
  return chords;
}
