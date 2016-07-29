#include <vector>
#include <iostream>

#include "Plotting/Plotting.hpp"
#include "SampleSorter/Tempo.hpp"
#include "SampleSorter/SpectralProcessing.hpp"

double Tempo::tempoToSeconds(double tempo) {
  return 1/tempo;
}

double Tempo::tempoToSamples(double tempo, long sampleRate) {
  return secondsToSamples(tempoToSeconds(tempo), sampleRate);
}

double Tempo::secondsToSamples(double seconds, long sampleRate) {
  return seconds * sampleRate;
}

double Tempo::samplesToSeconds(double samples, long sampleRate) {
  return samples/double(sampleRate);
}

double Tempo::tempoToBins(double tempo, long hopSize, long sampleRate) {
  return sampleRate/(tempo * hopSize);
}

double Tempo::binsToSamples(long bin, long hopSize) {
  return bin * hopSize;
}

double Tempo::binsToSeconds(long bin, long hopSize, long sampleRate) {
  return samplesToSeconds(binsToSamples(bin, hopSize), sampleRate);
}


double Tempo::tempoValue(double tempo, 
                         double theOne,
                         std::vector<double> onsets, 
                         long hopSize,
                         long sampleRate) {
  double value = 0;
  for (long i = 0; i < onsets.size(); i++) {
    double distanceFromBeat = (i - theOne)/tempoToBins(tempo, hopSize, sampleRate);
    distanceFromBeat = std::fmod(distanceFromBeat, 1);
    while (distanceFromBeat < 0)
      distanceFromBeat += 1;
    //if (distanceFromBeat > 0.5)
      //distanceFromBeat = 1 - distanceFromBeat;

    value += distanceFromBeat * onsets[i];
  }

  return value;//tempoToBins(tempo, hopSize, sampleRate);
}


double Tempo::correlationTempo(std::vector<double> onsets, 
                               long hopSize, 
                               long sampleRate) {
  //Plotting::plotVector(onsets, hopSize);

  // filter out frequencies less than 2 per clip
  // Clips with tempo must contain at least 2 beats
  double highPass = 2./binsToSeconds(onsets.size(), hopSize, sampleRate);
  // Also beat cannot be below 1/8 persecond
  highPass = std::max(highPass, 1/4.);

  // No tempos will exist above 1000bpm
  double lowPass = 1000/60.;

  std::vector<double> correlation = SpectralProcessing::autoCorrelation(onsets,
                                                    sampleRate/hopSize,
                                                    highPass,
                                                    lowPass);

  // only take first half of correlation
  long fftSize = correlation.size()/4 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(correlation.size()/2,
                                              correlation.data(),
                                              fft,
                                              FFTW_ESTIMATE);

  //Window the correlation
  for (long i = 0; i < correlation.size()/2; i++) {
    correlation[i] = correlation[i] * SpectralProcessing::hammingWindow(i, correlation.size()/2);
  }

  fftw_execute(fftPlan);
  fftw_destroy_plan(fftPlan);

  std::vector<std::pair<double, double> > peaks = 
    SpectralProcessing::findPeaks(fft, fftSize, sampleRate/hopSize);

  //Plotting::plotPair(peaks);

  fftw_free(fft);

  double max = 0;
  double tempo = 0;
  for (long i = 0; i < peaks.size(); i++) {
    if (peaks[i].second > max) {
      max = peaks[i].second;
      tempo = peaks[i].first;
    }
  }

  std::cout << "tempo = " << tempo * 60 << std::endl;
  return tempo;
}

std::pair<double, double> Tempo::fineTuneTempo(double tempo, 
    std::vector<double> onsets, 
    double BPSrange, 
    double BPSstepSize,
    long hopSize,
    long sampleRate) {

  std::vector<double> tempoValues;

  //Plotting::plotVector(onsets);

  double minTempoValue = onsets.size();
  double minTempo = tempo;
  double minOne = 0;
  for (double tempoAdjustment = -BPSrange; 
       tempoAdjustment <=BPSrange; 
       tempoAdjustment += BPSstepSize) {
    double trialTempo = tempo + tempoAdjustment;
    // adjust it by size of beats to bins -> less bins -> easier
    double thisMinValue = onsets.size();
    double thisMinTempo = 0;
    for (double theOne = 0; 
         theOne < tempoToBins(trialTempo, hopSize, sampleRate); 
         theOne ++) {
      double newTempoValue = tempoValue(trialTempo, theOne, onsets, hopSize, sampleRate);
      if (newTempoValue < minTempoValue) {
        minTempoValue = newTempoValue;
        minTempo = trialTempo;
        minOne = theOne;
      }

      if (newTempoValue < thisMinValue) {
        thisMinValue = newTempoValue;
        thisMinTempo = trialTempo;
      }
    }
    tempoValues.push_back(thisMinValue);
  }

  //Plotting::plotVector(tempoValues);
  std::cout << "tuned tempo = "  << minTempo*60 << std::endl;
  std::cout << "min one = " << binsToSeconds(minOne, hopSize, sampleRate) << std::endl;
  return std::make_pair(minTempo, binsToSeconds(minOne, hopSize, sampleRate));
}
