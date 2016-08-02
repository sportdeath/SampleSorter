#include <vector>
#include <iostream>

#include "SampleSorter/Sample.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/TimeDomainProcessing.hpp"
#include "SampleSorter/Tempo.hpp"
#include "Plotting/Plotting.hpp"
#include "SampleSorter/EqualLoudness.hpp"

Sample::Sample(std::string file_) {
  file = file_;
}

void Sample::process() {
  std::cout << "Processing '" << getName() << "'" << std::endl;
  tune();
  findBeat();
  findChords();
}


std::string Sample::getFile() {
  return file;
}

void Sample::tune() {
  std::vector<std::vector<double> > a = getWaves();
  EqualLoudness::filter(a, a, getSampleRate());
  TimeDomainProcessing::normalizeEnergy(a, a, getSampleRate());
  Octave o(a, 1200, getSampleRate(), tuningCents);
  tuningCents = o.tune();
  std::cout << "tuned by " << tuningCents << " cents" << std::endl;
  return;
}

void Sample::findBeat() {
  std::vector<std::vector<double> > a = getWaves();

  long hopSize = 1024;
  long windowRatio = 2;

  std::vector<double> onsets = SpectralProcessing::onsetEnergy(a, hopSize, windowRatio);

  tempo = Tempo::correlationTempo(onsets, hopSize, getSampleRate());
  std::pair<double, double> tempoOne = Tempo::fineTuneTempo(tempo, onsets, 10/60., 0.05/60., hopSize, getSampleRate());

  tempo = tempoOne.first;
  theOne = tempoOne.second;

  std::cout << "the one (with tuning): " << getTheOneWithTuning() << std::endl;
  std::cout << "Tempo (before tuning): " << 60*getRawBeat() << std::endl;
  std::cout << "Tempo (with tuning): " << 60*getBeatWithTuning() << std::endl;
}

void Sample::findChords() {
  long windowSize = Tempo::tempoToSamples(tempo, getSampleRate());
  long octaveSize = 12;

  std::vector<double> window(windowSize);

  // set up fourier transform
  long fftSize = windowSize/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
                                           window.data(),
                                           fft,
                                           FFTW_ESTIMATE);

  std::vector<std::vector<double> > a = getWaves();
  EqualLoudness::filter(a, a, getSampleRate());
  TimeDomainProcessing::normalizeEnergy(a, a, getSampleRate());

  long theOneSamples = Tempo::secondsToSamples(theOne, getSampleRate());
  long maxWindow = (a[0].size() - theOneSamples)/windowSize;

  std::vector<Octave> octaves(maxWindow);

  for (long hop = 0; hop < maxWindow; hop++) {
    octaves[hop] = Octave(octaveSize);
    for (long channel = 0; channel < a.size(); channel ++) {
      for (long i = 0; i < windowSize; i++) {
        window[i] = a[channel][theOneSamples + i + hop * windowSize];
        window[i] = window[i] * SpectralProcessing::hammingWindow(i, windowSize);
      }

      fftw_execute(fftPlan);
      Octave channelOctave(fft, fftSize, octaveSize, getSampleRate(), tuningCents);

      octaves[hop].add(octaves[hop], channelOctave);
    }
  }


  fftw_free(fft);
  fftw_destroy_plan(fftPlan);

  for (long i = 0; i < octaves.size(); i++) {
    std::cout << "@" << Tempo::binsToSeconds(i, windowSize, getSampleRate()) + theOne <<": " 
      << std::endl;
      //<< octaves[i].mostSimilarChord() << std::endl;
    octaves[i].plot();
  }
}

bool Sample::isHarmonic() {
  return isHarmonic_;
}

bool Sample::hasBeat() {
  return hasBeat_;
}

double Sample::getRawBeat() {
  return tempo;
}

double Sample::getBeatWithTuning() const {
  return tempo * pow(2., tuningCents/1200.);
}
double Sample::getTheOneWithTuning() const {
  return theOne / pow(2., tuningCents/1200.);
}

std::vector<Octave> Sample::getChords() {
  return chords;
}

bool Sample::isCompatible(const Sample & other) {
  // get larger tempo
  double larger, smaller;
  larger = std::max(getBeatWithTuning(), other.getBeatWithTuning());
  smaller = std::min(getBeatWithTuning(), other.getBeatWithTuning());
  long ratio = std::round(larger/smaller);
  double tuningSteps = 12 * std::log2(larger/(ratio * smaller));
  //double integerDeviation = std::fmod(tuningSteps, 1.);
  double integerDeviation = tuningSteps - std::round(tuningSteps);

  if (std::abs(integerDeviation) < 0.05) {
    std::cout << "'" << getName() << "' is compatible with '" << other.getName() <<"'" << std::endl;
    std::cout << "Tempos: " << 60*getBeatWithTuning() <<", " << 60*other.getBeatWithTuning() << std::endl;
    std::cout << "Tuning cents: " << tuningCents <<", " << other.tuningCents << std::endl;
    std::cout << "Tuning ratio: " << std::round(tuningSteps) << " steps and " << integerDeviation * 100 << " cents" << std::endl << std::endl;
  }

  return integerDeviation < 0.05;
}

  // find whether that ratio is an even number of cents

//bool Sample::isSimilar(const & Sample other) {
  //// coarse tune so their tempos are equal
  //// or are integer ratios of each other
  //// get similarity of all octaves.
  //// Look for large chains
//}
