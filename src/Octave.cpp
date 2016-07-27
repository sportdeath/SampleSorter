#include <vector>
#include <cmath>
#include <complex>

#include <fftw3.h>

#include "gnuplot-iostream.h"

#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/Octave.hpp"

Octave::Octave(std::vector< std::vector<double> > audio, 
       long bins,
       long sampleRate) {

  BASE_OFFSET = 1. - fmod(log2(BASE_FREQ), 1);
  CENTS_PER_BIN = CENTS_PER_OCTAVE/double(bins);
  BINS_PER_SEMITONE = bins/double(SEMITONES_PER_OCTAVE);

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

    // for all peaks, get peak frequency and add peak
    std::vector< std::pair<double, double> > peaks = 
      SpectralProcessing::findPeaks(fft, fftSize, sampleRate);

    // Add all the peaks to the spectrogram
    for (long i = 0; i < peaks.size(); i++) {
      addPeak(peaks[i].first, peaks[i].second);
    }
  }

  fftw_destroy_plan(fftPlan);
  fftw_free(fft);

  // Normalize so that magnitude is 1
  double magnitude = 0;
  for (long i = 0; i < spectrogram.size(); i++) {
    magnitude += spectrogram[i] * spectrogram[i];
  }
  magnitude = sqrt(magnitude);
  for (long i = 0; i < spectrogram.size(); i++) {
    spectrogram[i] = spectrogram[i]/magnitude;
  }
}

void Octave::addPeak(double peakFreq, double peakValue) {
  double freqMod = fmod(log2(peakFreq) + BASE_OFFSET, 1);
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
  std::vector<std::pair<double, double> > xyPoints;
  for (long i = 0; i < spectrogram.size(); i++) {
    xyPoints.push_back(std::make_pair(i, spectrogram[i]));
  }
  Gnuplot gp;
  gp << "plot" << gp.file1d(xyPoints) << "w l" << std::endl;
  std::cin.get();
}
