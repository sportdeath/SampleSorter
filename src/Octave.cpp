#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/Octave.hpp"
#include "Plotting/Plotting.hpp"

const std::vector<double> Octave::harmonicSeries( {0, 5, 4, 3, 2, 1, 6, 1, 2, 3, 4, 5});
const std::vector<std::string> Octave::noteNames( {"A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"});

const Octave Octave::unison (std::vector<double>({1,0,0,0,0,0,0,0,0,0,0,0}));
const Octave Octave::minorSecond (std::vector<double>({1,1,0,0,0,0,0,0,0,0,0,0}));
const Octave Octave::majorSecond (std::vector<double>({1,0,1,0,0,0,0,0,0,0,0,0}));
const Octave Octave::minorThird  (std::vector<double>({1,0,0,1,0,0,0,0,0,0,0,0}));
const Octave Octave::majorThird  (std::vector<double>({1,0,0,0,1,0,0,0,0,0,0,0}));
const Octave Octave::perfectFifth(std::vector<double>({1,0,0,0,0,0,0,1,0,0,0,0}));

const Octave Octave::majorTriad     (std::vector<double>({1,0,0,0,1,0,0,1,0,0,0,0}));
const Octave Octave::augmentedTriad (std::vector<double>({1,0,0,0,1,0,0,0,1,0,0,0}));
const Octave Octave::minorTriad     (std::vector<double>({1,0,0,1,0,0,0,1,0,0,0,0}));
const Octave Octave::diminishedTriad(std::vector<double>({1,0,0,1,0,0,1,0,0,0,0,0}));

const std::vector<Octave> Octave::chords(std::vector<Octave>({unison,
                                                              minorSecond,
                                                              majorSecond,
                                                              minorThird,
                                                              majorThird,
                                                              perfectFifth,
                                                              majorTriad,
                                                              augmentedTriad,
                                                              diminishedTriad
                                                              }));

Octave::Octave(long bins) 
{
  CENTS_PER_BIN = CENTS_PER_OCTAVE/double(bins);
  BINS_PER_SEMITONE = bins/double(SEMITONES_PER_OCTAVE);

  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);
}

Octave::Octave(std::vector<double> spec) : Octave(spec.size()) {
  spectrogram = spec;
  normalize();
}

Octave::Octave(std::vector< std::vector<double> > audio, 
       long bins,
       long sampleRate, 
       long tuningCents) : Octave(bins) {

  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);

  long fftSize = audio[0].size()/2 + 1;

  std::vector<double> windows(audio[0].size());
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windows.size(),
                                 windows.data(),
                                 fft,
                                 FFTW_ESTIMATE |
                                 FFTW_DESTROY_INPUT);

  for (long channel = 0; channel < audio.size(); channel ++) {
    // Window the audio
    SpectralProcessing::hammingWindow(windows, audio[channel]);

    // Take the fourier transform
    fftw_execute(fftPlan);
    
    // Make octave with fourier transform
    Octave a(fft, fftSize, bins, sampleRate, tuningCents);
    this -> add(a);
  }

  fftw_destroy_plan(fftPlan);
  fftw_free(fft);
}

Octave::Octave(fftw_complex * fft,
       long fftSize,
       long bins,
       long sampleRate, 
       long tuningCents) : Octave(bins) {

  double baseOffset = 1. - fmod(log2(BASE_FREQ) - tuningCents/1200., 1);

  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);

  // for all peaks, get peak frequency and add peak
  std::vector< std::pair<double, double> > peaks = 
    SpectralProcessing::findPeaks(fft, fftSize, sampleRate);

  // Add all the peaks to the spectrogram
  for (long i = 0; i < peaks.size(); i++) {
    addPeak(peaks[i].first, peaks[i].second, baseOffset);
  }

  // Normalize so that magnitude is 1
  normalize();
}

void Octave::normalize() {
  double magnitude = 0;
  for (long i = 0; i < spectrogram.size(); i++) {
    magnitude += spectrogram[i] * spectrogram[i];
  }
  magnitude = sqrt(magnitude);
  for (long i = 0; i < spectrogram.size(); i++) {
    spectrogram[i] = spectrogram[i]/magnitude;
  }
}


void Octave::addPeak(double peakFreq, double peakValue, double baseOffset) {
  if (peakFreq > 5000) {
    // Humans do not consider pitches
    // greater than 5000 to have harmonic content
    return;
  }

  double freqMod = fmod(log2(peakFreq) + baseOffset, 1);
  while (freqMod < 0) {
    freqMod += 1;
  } 

  double desiredBin = freqMod * spectrogram.size();

  int leftBin = floor(desiredBin);
  int rightBin = (leftBin + 1) % spectrogram.size();

  double rightPercent = desiredBin - leftBin;
  double leftPercent = 1 - rightPercent;

  double rightAmp = rightPercent * peakValue;
  double leftAmp = leftPercent * peakValue;

  spectrogram[leftBin] += leftAmp;
  spectrogram[rightBin] += rightAmp;
}

void Octave::rotate(long bins) {
  while (bins < 0) {
    bins += spectrogram.size();
  }

  std::rotate(spectrogram.begin(), spectrogram.end() - bins, spectrogram.end());
}

double Octave::tuningValue() {
  // Value at bin * distance of bin from 
  double output = 0;
  for (long bin = 0; bin < spectrogram.size(); bin ++) {
    double distanceFromSemitone = std::fmod(bin/(spectrogram.size()/double(SEMITONES_PER_OCTAVE)), 1.);
    if (distanceFromSemitone > 0.5) {
      distanceFromSemitone = 1 - distanceFromSemitone;
    }
    output += distanceFromSemitone * spectrogram[bin];
  }
  double sum = 0;
  return output;
}

long Octave::tune() {
  // rotate by -50 cents
  normalize();
  double max = std::ceil(BINS_PER_SEMITONE/2.);
  rotate(-max);
  double minTuningValue = spectrogram.size();
  long minBin = -max;

  // rotate by +1 cent 100x
  for (long bin = -max; bin <= max; bin++) {
    double currentTuningValue = tuningValue();
    if (currentTuningValue < minTuningValue) {
      minTuningValue = currentTuningValue;
      minBin = bin;
    }
    rotate(1);
  }
  // return to minimum and return shift
  rotate(minBin-(max+1));
  return minBin * CENTS_PER_BIN;
}

void Octave::plot() {
  Plotting::plotVector(spectrogram, 1/BINS_PER_SEMITONE);
}

double Octave::difference(const Octave & other) {
  if (spectrogram.size() != other.spectrogram.size()) {
    std::cout << spectrogram.size() <<"!="<< other.spectrogram.size() << std::endl;
    throw 1;
  }

  double output = 0;
  for (long i = 0; i < spectrogram.size(); i++) {
    output += spectrogram[i] * other.spectrogram[i];
  }

  output = std::max(std::min(output, 1.), 0.);

  output = std::acos(output) * 2/M_PI;

  return output;
}

double Octave::similarity(const Octave & other) {
  return 1 - difference(other);
}

void Octave::add(const Octave & other) {
  if (spectrogram.size() != other.spectrogram.size()) {
    throw 1;
  }


  for (long i = 0; i < spectrogram.size(); i++) {
    spectrogram[i] += other.spectrogram[i];
  }

  normalize();
}

Octave Octave::mean(
    std::vector<Octave>::iterator begin,
    std::vector<Octave>::iterator end) {
  std::vector<double> totalSpectrogram(begin -> spectrogram.size());
  std::fill(totalSpectrogram.begin(), totalSpectrogram.end(), 0);
  for (auto it = begin; it != end; it++) {
    for (long i = 0; i < it -> spectrogram.size(); i++) {
      totalSpectrogram[i] += it -> spectrogram[i];
    }
  }
  return Octave(totalSpectrogram);
}

double Octave::selfSimilarity(
    std::vector<Octave>::iterator begin,
    std::vector<Octave>::iterator end) {
  const Octave m = mean(begin,end);
  double totalSimilarity = 0;
  for (auto it = begin; it != end; it++) {
    totalSimilarity += (*it).similarity(m);
  }
  return totalSimilarity;
    ///double(std::distance(begin,end));
}

double Octave::dissonance() {
  double total = 0;
  for (long root = 0; root < spectrogram.size(); root++) {
    double subTotal = 0;
    for (long harmNum = 1; harmNum < spectrogram.size(); harmNum++) {
      long index = (root + harmNum) % spectrogram.size();
      subTotal += spectrogram[index] * harmonicSeries[harmNum];
    }
    total += spectrogram[root] * subTotal;
  }
  return total;
}
    
double Octave::dissonance(
    std::vector<Octave>::iterator begin,
    std::vector<Octave>::iterator end) {
  Octave m = mean(begin, end);
  return m.dissonance();
}

std::string Octave::mostSimilarChord() {
  double maxSimilarity = 0;
  std::string best;
  for (long rotation = 0; rotation < 12; rotation ++) {
    for (long chord = 0; chord < chords.size(); chord++) {
      double newSimilarity = similarity(chords[chord]);
      if (newSimilarity > maxSimilarity) {
        maxSimilarity = newSimilarity;
        best = noteNames[rotation] + " " + std::to_string(chord);
      }
    }
    rotate(1);
  }
  return best;
}
