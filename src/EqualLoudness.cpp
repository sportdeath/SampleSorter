#include <vector>
#include <cmath>
#include <iostream>
#include <ciso646>

#include <fftw3.h>

#include "SampleSorter/EqualLoudness.hpp"

double EqualLoudness::AWeightingAmp(double freq) {
  if (freq > 30000 or freq < 10)
    return 0;
  double freqSquared = freq*freq;
  double top = 12194*12194*freqSquared*freqSquared;
  double bot1 = 20.6*20.6 + freqSquared;
  double bot2 = 107.7*107.7 + freqSquared;
  double bot3 = 737.9*737.9 + freqSquared;
  double bot4 = 12194*12194 + freqSquared;
  return top/(bot1*std::sqrt(bot2*bot3)*bot4);
}

double EqualLoudness::CWeightingAmp(double freq) {
  if (freq > 30000 or freq < 10)
    return 0;
  double freqSquared = freq*freq;
  double top = 148840000*freqSquared;
  double bot1 = 424.36 + freqSquared;
  double bot2 = 1488840000 + freqSquared;
  return top/(bot1*bot2);
}

void EqualLoudness::filter(
    std::vector<std::vector<double> > & output,
    std::vector<std::vector<double> > & audio,
    long sampleRate) {
  long fftSize = audio[0].size()/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

  output.resize(audio.size());

  for (long channel = 0; channel < audio.size(); channel ++) {
    fftw_plan fftPlanFor = fftw_plan_dft_r2c_1d(audio[channel].size(),
                                             audio[channel].data(),
                                             fft,
                                             FFTW_ESTIMATE);
    // take fourier transform
    fftw_execute(fftPlanFor);
    for (long bin = 0; bin < fftSize; bin++) {
      // get freq of bin
      double freq = bin * sampleRate/audio[channel].size();
      // multiply by weighting
      fft[bin][0] *= AWeightingAmp(freq)/audio[channel].size();
      fft[bin][1] *= AWeightingAmp(freq)/audio[channel].size();
      //fft[bin][0] *= 1./audio[channel].size();
      //fft[bin][1] *= 1./audio[channel].size();
    }

    fftw_destroy_plan(fftPlanFor);

    output[channel].resize(audio[channel].size());
    fftw_plan fftPlanInv = fftw_plan_dft_c2r_1d(audio[channel].size(),
                                             fft,
                                             output[channel].data(),
                                             FFTW_ESTIMATE);

    //inverse
    fftw_execute(fftPlanInv);
    fftw_destroy_plan(fftPlanInv);
  }
  fftw_free(fft);
}
