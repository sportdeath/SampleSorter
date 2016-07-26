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
  CENTS_PER_BIN = CENTS_PER_OCTAVE/bins;

  spectrogram.resize(bins);
  std::fill(spectrogram.begin(), spectrogram.end(), 0);


  for (long channel = 0; channel < audio.size(); channel ++) {
    // Window the audio
    std::vector<double> windows = SpectralProcessing::hammingWindow(audio[channel]);

    long fftSize = windows.size()/2 + 1;
    fftw_complex fft[fftSize];
    fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windows.size(),
                                   windows.data(),
                                   fft,
                                   FFTW_DESTROY_INPUT);

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

  std::cout << "rotating by " << bins << std::endl;

  // shift vector by floor[cents/CENTS_PER_BIN]
  //long integerShift = cents/CENTS_PER_BIN;
  std::rotate(spectrogram.begin(), spectrogram.end() - bins, spectrogram.end());

  //double percentShift = std::fmod(cents/double(CENTS_PER_BIN), 1.);
  //std::cout << percentShift << std::endl;

  //double lastEntry = spectrogram[spectrogram.size() - 1];

  //for (long i = spectrogram.size() - 1; i > 0; i--) {
    //spectrogram[i] = spectrogram[i] * (1. - percentShift) + spectrogram[i - 1] * (percentShift);
  //}
  //spectrogram[0] = spectrogram[0] * (1. - percentShift) + lastEntry * (percentShift);
  

}

double Octave::tuningValue() {
  // Value at bin * distance of bin from 
  double output = 0;
  for (long bin = 0; bin < spectrogram.size(); bin ++) {
    double distanceFromSemitone = std::fmod(bin/(spectrogram.size()/double(SEMITONES_PER_OCTAVE)), 1.);
    if (distanceFromSemitone > 0.5) {
      distanceFromSemitone = 1 - distanceFromSemitone;
    }
    std::cout << "bin: " << bin << " distanceFromSemitone: " << distanceFromSemitone << std::endl;
    output += distanceFromSemitone * spectrogram[bin];
  }
  std::cout << "tuningValue: " << output << std::endl;
  plot();
  double sum = 0;
  return output;
}

long Octave::tune() {
  // rotate by -50 cents
  // rotate by +1 cent 100x
  // return to minimum and return shift
  double max = std::ceil(spectrogram.size()/double(SEMITONES_PER_OCTAVE)/2.);
  rotate(-max);
  double minTuningValue = spectrogram.size();
  long minBin = 0;
  for (long bin = -max; bin <= max; bin++) {
    rotate(1);
    double currentTuningValue = tuningValue();
    if (currentTuningValue < minTuningValue) {
      minTuningValue = currentTuningValue;
      minBin = bin;
    }
  }
  //rotate(-CENTS_PER_SEMITONE/2);
  //rotate(minCent);
  return minBin;
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
