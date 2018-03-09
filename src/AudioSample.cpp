#include <vector>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/TimeDomainProcessing.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/EqualLoudness.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/Units.hpp"
#include "SampleSorter/Tempo.hpp"

#include "Plotting/Plotting.hpp"

// 20% error max
const double AudioSample::TEMPO_PERCENTAGE_ERROR = 0.2;
const double AudioSample::TEMPO_STEPS = 2000;
const double AudioSample::TEMPO_ONE_STEPS = 1000;

AudioSample::AudioSample() {
  tuningCents = 0;
}

AudioSample::AudioSample(
    std::vector<std::vector<double> > & audio,
    long _sampleRate
    ) : sampleRate(_sampleRate) {

  totalSeconds = Units::samplesToSeconds(audio[0].size(), sampleRate);
  // filter the audio
  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);

  tune(filteredAudio);
  findFundemental(filteredAudio);
  // do not use the filtered audio here
  // the information is useful
  findBeat(audio);

  // no more chords
  // findChords(filteredAudio);
  // chords.resize(1);
  // chords[0] = octave;
}

double AudioSample::getTotalSeconds() const {
  return totalSeconds;
}

AudioSample::AudioSample(
    long tuningCents_,
    short fundemental_,
    double rawBeat,
    double theOneBin,
    double totalSeconds_,
    long sampleRate,
    std::vector<Octave> chords_,
    Octave octave_) :
  tempo(rawBeat, theOneBin, sampleRate),
  sampleRate(sampleRate)
{
  totalSeconds = totalSeconds_;
  tuningCents = tuningCents_;
  fundemental = fundemental_;
  chords = chords_;
  octave = octave_;
}

void AudioSample::tune(std::vector<std::vector<double> > & audio) {
  Octave oct(audio, 1200., sampleRate);
  tuningCents = oct.tune();
}

void AudioSample::findFundemental(std::vector<std::vector<double> > & audio) {
  octave = Octave(audio, 12., sampleRate, tuningCents);
  fundemental = octave.getMax();
}

void AudioSample::findBeat(std::vector<std::vector<double> > & audio) {
  tempo = Tempo(
      audio, 
      sampleRate, 
      TEMPO_PERCENTAGE_ERROR,
      TEMPO_STEPS,
      TEMPO_ONE_STEPS
      );
}

void AudioSample::findChords(std::vector<std::vector<double> > & audio) {
  TimeDomainProcessing::unitEnergyPerBeat(
      audio, audio, tempo.getTempo(), sampleRate);

  int windowSize = Units::tempoToSamples(tempo.getTempo(), sampleRate);
  long startSample = Units::secondsToSamples(tempo.getTheOne(), sampleRate);

  int numChords = 
    (Units::secondsToSamples(totalSeconds, sampleRate) - startSample)/windowSize;

  chords.resize(numChords);

  std::vector<double> window(windowSize);

  // set up fft
  long fftSize = windowSize/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
                                           window.data(),
                                           fft,
                                           FFTW_ESTIMATE);

  for (int i = 0; i < numChords; i++) {

    Octave chord;

    // for each channel
    for (long channel = 0; channel < audio.size(); channel ++) {
      // window this section
      for (long i = 0; i < windowSize; i++) {
        window[i] = audio[channel][i + startSample];
        window[i] = window[i] * SpectralProcessing::hammingWindow(i, windowSize);
      }

      // take fourier transform of this section
      fftw_execute(fftPlan);
      // turn it into a chord
      Octave channelChord(fft, fftSize, 12, sampleRate, tuningCents);

      // add to channel
      chord.add(chord, channelChord);
    }

    // add to chords
    chords[i] = chord;
  }

  fftw_free(fft);
  fftw_destroy_plan(fftPlan);
}

long AudioSample::getTuningCents() const {
  return tuningCents;
}
short AudioSample::getFundemental() const {
  return fundemental;
}
double AudioSample::getTuningCentsFreqRatio() const {
  return std::pow(2., tuningCents/1200.);
}
double AudioSample::getBeatRaw() const {
  return tempo.getTempo();
}
double AudioSample::getBeatWithTuning() const {
  return tempo.getTempo() * getTuningCentsFreqRatio();
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
Octave AudioSample::getOctave() const {
  return octave;
}
long AudioSample::getSampleRate() const {
  return sampleRate;
}
double AudioSample::getLastBeatSeconds() const {
  int numBeats = 
    std::floor((totalSeconds - getTheOneRaw()) * getBeatRaw());
  return (numBeats / getBeatRaw()) + getTheOneRaw();
}
double AudioSample::energy(std::vector<std::vector<double> > & audio) const {
    double energy = 0;
    for (long i = 0; i < audio[0].size(); i++) {
        double sample = 0;
        for (long channel = 0; channel < audio.size(); channel++) {
            sample += audio[channel][i];
        }
        energy += sample * sample;
    }
    return energy;
}
