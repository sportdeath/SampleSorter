#include <vector>
#include <iostream>

#include "SampleSorter/Sample.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/Tempo.hpp"
#include "Plotting/Plotting.hpp"
#include "SampleSorter/EqualLoudness.hpp"

Sample::Sample(std::string file_) {
  file = file_;
}

void Sample::process() {
  tune();
  findBeat();
  findChords();
}


std::string Sample::getFile() {
  return file;
}

void Sample::tune() {
  std::vector<std::vector<double> > a = getWaves();
  a = EqualLoudness::filter(a, getSampleRate());
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
}

void Sample::findChords() {
  long hopSize = 1024;
  long octaveSize = 12;

  std::vector<double> window(hopSize);

  // set up fourier transform
  long fftSize = hopSize/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(hopSize,
                                           window.data(),
                                           fft,
                                           FFTW_MEASURE);

  std::vector<std::vector<double> > a = getWaves();
  long maxHop = a[0].size()/hopSize;

  std::vector<Octave> octaves(maxHop);

  for (long hop = 0; hop < maxHop; hop++) {
    octaves[hop] = Octave(octaveSize);
    for (long channel = 0; channel < a.size(); channel ++) {
      for (long i = 0; i < hopSize; i++) {
        window[i] = a[channel][i + hop * hopSize];
        window[i] = window[i] * SpectralProcessing::hammingWindow(i, hopSize);
      }

      fftw_execute(fftPlan);
      Octave channelOctave(fft, fftSize, octaveSize, getSampleRate(), tuningCents);

      octaves[hop].add(channelOctave);
    }
  }


  fftw_free(fft);
  fftw_destroy_plan(fftPlan);

  // could use hopsize so that its a nice multiple...
  long searchSize = getSampleRate()/double(getRawBeat() * hopSize);

  std::vector<double> dissonances(searchSize);

  for (long searchPos = 0; searchPos < searchSize; searchPos ++) {
    std::vector<Octave>::iterator begin = octaves.begin();
    std::vector<Octave>::iterator end = begin + searchPos;
    double dissonance = 0;
    while (end < octaves.end()) {
      dissonance += Octave::dissonance(begin, end);
      begin = end;
      end = begin + searchSize;
    }
    end = octaves.end();
    dissonance += Octave::dissonance(begin, end);
    dissonances[searchPos] = dissonance;
  }

  double minimum = dissonances[0];
  long minIndex =0;
  for (long i = 0; i < searchSize; i++) {
    if (dissonances[i] < minimum) {
      minimum = dissonances[i];
      minIndex = i;
    }
  }

  std::cout << "The 1 found at " << (minIndex * hopSize)/double(getSampleRate()) << " seconds" << std::endl;

  //Plotting::plotVector(dissonances);

  std::vector<Octave>::iterator begin = octaves.begin();
  std::vector<Octave>::iterator end = begin + minIndex;
  while (end < octaves.end()) {
    Octave m = Octave::mean(begin, end);
    chords.push_back(m);
    begin = end;
    end = begin + searchSize;
  }
  end = octaves.end();
  Octave m = Octave::mean(begin, end);
  chords.push_back(m);
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

double Sample::getBeatWithTuning() {
  return tempo;
}

std::vector<Octave> Sample::getChords() {
  return chords;
}

//bool Sample::isSimilar(const & Sample other) {
  //// coarse tune so their tempos are equal
  //// or are integer ratios of each other
  //// get similarity of all octaves.
  //// Look for large chains
//}
