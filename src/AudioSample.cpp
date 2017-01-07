#include <vector>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/TimeDomainProcessing.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/EqualLoudness.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/Units.hpp"
#include "SampleSorter/TempoFunction.hpp"

AudioSample::AudioSample() {
  tuningCents = 0;
}

AudioSample::AudioSample(
    std::vector<std::vector<double> > & audio, 
    long sampleRate) {

  tune(audio, sampleRate);
  findBeat(audio, sampleRate);
  findChords(audio, sampleRate);
}

double AudioSample::getTotalSeconds() const {
  return tempo.getTotalSeconds();
}

AudioSample::AudioSample(
    long tuningCents_,
    double rawBeat,
    double theOne,
    double totalSeconds,
    std::vector<Octave> chords_) {
  tuningCents = tuningCents_;
  tempo = TempoFunction(rawBeat, theOne, totalSeconds);
  chords = chords_;
}

void AudioSample::tune(std::vector<std::vector<double> > & audio, long sampleRate) {
  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);

  Octave oct(filteredAudio, 1200., sampleRate);
  tuningCents = oct.tune();
}

void AudioSample::findBeat(std::vector<std::vector<double> > & audio, long sampleRate) {
  double percentageError = 0.5; // 10%
  double percentageStep = 0.0001; // 0.1%

  // tempo drift has resolution up to 4 seconds
  size_t degrees = Units::samplesToSeconds(audio[0].size(), sampleRate)/4.;

  tempo = TempoFunction(audio, 
                        sampleRate, 
                        degrees, 
                        percentageError, 
                        percentageStep);
}

void AudioSample::findChords(std::vector<std::vector<double> > & audio, long sampleRate) {
  // filter the audio
  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);
  TimeDomainProcessing::unitEnergyPerBeat(
      filteredAudio, filteredAudio, tempo.getAvgTempo(), sampleRate);

  long k = 0;
  double startBeat = 0;
  double endBeat = std::min(tempo.getTotalSeconds(), tempo.getKthBeat(k));

  chords.resize(0);

  do {
    // set up window
    long startSample = Units::secondsToSamples(startBeat, sampleRate);
    long endSample = Units::secondsToSamples(endBeat, sampleRate);
    long windowSize = endSample - startSample;
    if (windowSize > 0 and startSample >= 0) {
      std::vector<double> window(windowSize);

      // set up fft
      long fftSize = windowSize/2 + 1;
      fftw_complex * fft;
      fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
      fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
                                               window.data(),
                                               fft,
                                               FFTW_ESTIMATE);

      Octave chord;

      // for each channel
      for (long channel = 0; channel < filteredAudio.size(); channel ++) {
        // window this section
        //std::cout << "before" << std::endl;
        for (long i = 0; i < windowSize; i++) {
          window[i] = filteredAudio[channel][i + startSample];
          window[i] = window[i] * SpectralProcessing::hammingWindow(i, windowSize);
        }
        //std::cout << "after" << std::endl;

        // take fourier transform of this section
        fftw_execute(fftPlan);
        // turn it into a chord
        Octave channelChord(fft, fftSize, 12, sampleRate, tuningCents);

        // add to channel
        chord.add(chord, channelChord);
      }

      fftw_free(fft);
      fftw_destroy_plan(fftPlan);

      // add to chords
      chords.push_back(chord);
    }

    k += 1;
    startBeat = endBeat;
    endBeat = std::min(tempo.getTotalSeconds(), tempo.getKthBeat(k));
  } while (startBeat != tempo.getTotalSeconds());
}

long AudioSample::getTuningCents() const {
  return tuningCents;
}
double AudioSample::getTuningCentsFreqRatio() const {
  return std::pow(2., tuningCents/1200.);
}
double AudioSample::getBeatRaw() const {
  return tempo.getAvgTempo();
}
double AudioSample::getBeatWithTuning() const {
  return tempo.getAvgTempo() * getTuningCentsFreqRatio();
}
double AudioSample::getTheOneRaw() const {
  return tempo.getTheOne();
}
double AudioSample::getTheOneWithTuning() const {
  return tempo.getTheOne()/getTuningCentsFreqRatio();
}
std::vector<Octave> AudioSample::getChords() const {
  return chords;
}
