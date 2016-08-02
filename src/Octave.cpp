#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/MusicTheory.hpp"
#include "Plotting/Plotting.hpp"

Octave::Octave(long bins) {
  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);
}

Octave::Octave(std::vector<double> spec) {
  spectrogram = spec;
}

Octave::Octave(std::vector< std::vector<double> > audio, 
       long bins,
       long sampleRate, 
       long tuningCents) : Octave(bins) {

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
    //for (long i = 0; i < audio[channel].size(); i++) {
      //windows[i] = audio[channel][i];
    //}

    // Take the fourier transform
    fftw_execute(fftPlan);
    
    // Make octave with fourier transform
    Octave channelOctave(fft, fftSize, bins, sampleRate, tuningCents);
    add(*this, channelOctave);
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

  // for all peaks, get peak frequency and add peak
  std::vector< std::pair<double, double> > peaks = 
    SpectralProcessing::findPeaks(fft, fftSize, sampleRate);

  // Add all the peaks to the spectrogram
  for (long i = 0; i < peaks.size(); i++)
    addPeak(peaks[i].first, 
            peaks[i].second,
            baseOffset);
}

long Octave::getBins() const {
  return spectrogram.size();
}

double Octave::getBinsPerSemitone() const {
  return getBins()/double(SEMITONES_PER_OCTAVE);
}

double Octave::getCentsPerBin() const {
  return CENTS_PER_OCTAVE/double(getBins());
}

void Octave::addPeak(double peakFreq, double peakValue, double baseOffset) {
  if (peakFreq > 3000)
    // Humans do not consider pitches
    // greater than 5000 to have harmonic content
    return;

  double freqMod = fmod(log2(peakFreq) + baseOffset, 1);
  while (freqMod < 0) {
    freqMod += 1;
  } 

  // We want the spectral energy:

  double desiredBin = freqMod * getBins();

  long leftBin = floor(desiredBin);
  long rightBin = (leftBin + 1) % getBins();

  double rightPercent = desiredBin - leftBin;
  double leftPercent = 1 - rightPercent;

  double rightAmp = rightPercent * peakValue;
  double leftAmp = leftPercent * peakValue;

  spectrogram[leftBin] += leftAmp;
  spectrogram[rightBin] += rightAmp;
}

void Octave::rotate(long bins) {
  while (bins < 0) {
    bins += getBins();
  }

  std::rotate(spectrogram.begin(), spectrogram.end() - bins, spectrogram.end());
}

double Octave::tuningValue() const {
  // Value at bin * distance of bin from 
  double output = 0;
  for (long bin = 0; bin < getBins(); bin ++) {
    double distanceFromSemitone = bin/getBinsPerSemitone();
    distanceFromSemitone = std::abs(distanceFromSemitone 
                                    - std::round(distanceFromSemitone));

    output += distanceFromSemitone * spectrogram[bin];
  }
  return output;
}

long Octave::tune() {
  // rotate by -50 cents
  double max = std::ceil(getBinsPerSemitone()/2.);
  rotate(-max);

  double minTuningValue = tuningValue();
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
  return minBin * getCentsPerBin();
}

void Octave::plot() const {
  Plotting::plotVector(spectrogram, 1/getBinsPerSemitone());
}

void Octave::add(Octave & output, const Octave & other) const {
  if (getBins() != other.getBins() or getBins() != output.getBins()) {
    throw 1;
  }

  for (long bin = 0; bin < getBins(); bin++) {
    output.spectrogram[bin] = spectrogram[bin] + other.spectrogram[bin];
    output.spectrogram[bin] /= 2.;
  }
}
