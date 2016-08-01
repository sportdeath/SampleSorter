#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/MusicTheory.hpp"
#include "Plotting/Plotting.hpp"

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
       long tuningCents,
       bool quantize) : Octave(bins) {

  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);

  long fftSize = audio[0].size()/2 + 1;

  std::cout << fftSize << std::endl;

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
       long tuningCents,
       bool quantize) : Octave(bins) {

  double baseOffset = 1. - fmod(log2(BASE_FREQ) - tuningCents/1200., 1);

  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);

  // for all peaks, get peak frequency and add peak
  std::vector< std::pair<double, double> > peaks = 
    SpectralProcessing::findPeaks(fft, fftSize, sampleRate);

  double maxPeak = 0;
  for (long i = 0; i < peaks.size(); i++)
    maxPeak = std::max(maxPeak, peaks[i].second);

  // Add all the peaks to the spectrogram
  for (long i = 0; i < peaks.size(); i++)
    addPeak(peaks[i].first, peaks[i].second/double(fftSize), baseOffset, maxPeak/double(fftSize), quantize);

  Plotting::plotPair(peaks);
  plot();

  // Normalize so that magnitude is 1
  normalize();
}

void Octave::normalize() {
  double newMagnitude = 0;
  for (long i = 0; i < spectrogram.size(); i++) {
    newMagnitude += spectrogram[i] * spectrogram[i];
  }
  newMagnitude = sqrt(newMagnitude);
  for (long i = 0; i < spectrogram.size(); i++) {
    spectrogram[i] = spectrogram[i]/newMagnitude;
  }
  magnitude = newMagnitude * magnitude;
}


void Octave::addPeak(double peakFreq, double peakValue, double baseOffset, double maxPeak, bool quantize) {
  if (peakFreq > 1500)
    // Humans do not consider pitches
    // greater than 5000 to have harmonic content
    return;
  if (peakValue < 0.5*maxPeak)
    return;

  double freqMod = fmod(log2(peakFreq) + baseOffset, 1);
  while (freqMod < 0) {
    freqMod += 1;
  } 

  double desiredBin = freqMod * spectrogram.size();

  if (quantize) {
    long bin = std::fmod(std::round(desiredBin), spectrogram.size());

    spectrogram[bin] += peakValue;
  } else {
    long leftBin = floor(desiredBin);
    long rightBin = (leftBin + 1) % spectrogram.size();

    double rightPercent = desiredBin - leftBin;
    double leftPercent = 1 - rightPercent;

    double rightAmp = rightPercent * peakValue;
    double leftAmp = leftPercent * peakValue;

    spectrogram[leftBin] += leftAmp;
    spectrogram[rightBin] += rightAmp;
  }
}

void Octave::rotate(long bins) {
  while (bins < 0) {
    bins += spectrogram.size();
  }

  std::rotate(spectrogram.begin(), spectrogram.end() - bins, spectrogram.end());
}

double Octave::tuningValue() const {
  // Value at bin * distance of bin from 
  double output = 0;
  for (long bin = 0; bin < spectrogram.size(); bin ++) {
    double distanceFromSemitone = bin/(spectrogram.size()/double(SEMITONES_PER_OCTAVE));
    distanceFromSemitone = std::abs(distanceFromSemitone 
                                    - std::round(distanceFromSemitone));

    output += distanceFromSemitone * spectrogram[bin];
  }
  double sum = 0;
  return output;
}

long Octave::tune() {
  // rotate by -50 cents
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

void Octave::plot() const {
  Plotting::plotVector(spectrogram, 1/BINS_PER_SEMITONE);
}

double Octave::difference(const Octave & other) const {
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

double Octave::similarity(const Octave & other) const {
  return 1 - difference(other);
}

void Octave::add(const Octave & other) {
  if (spectrogram.size() != other.spectrogram.size()) {
    throw 1;
  }

  for (long i = 0; i < spectrogram.size(); i++) {
    spectrogram[i] = spectrogram[i] * magnitude 
                    + other.spectrogram[i] * other.magnitude;
  }
  
  magnitude = 1;

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
  //double total = 0;
  //for (long root = 0; root < spectrogram.size(); root++) {
    //double subTotal = 0;
    //for (long harmNum = 1; harmNum < spectrogram.size(); harmNum++) {
      //long index = (root + harmNum) % spectrogram.size();
      //subTotal += spectrogram[index] * harmonicSeries[harmNum];
    //}
    //total += spectrogram[root] * subTotal;
  //}
  //return total;
  return 0;
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
    for (long chord = 0; chord < MusicTheory::chords.size(); chord++) {
      double newSimilarity = similarity(MusicTheory::chords[chord].second);
      if (newSimilarity > maxSimilarity) {
        maxSimilarity = newSimilarity;
        best = MusicTheory::noteNames[rotation] + " " + MusicTheory::chords[chord].first;
      }
    }
    rotate(-1);
  }
  return best;
}
