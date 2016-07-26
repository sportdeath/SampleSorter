#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

#include "SampleSorter/SpectralProcessing.hpp"

double SpectralProcessing::hammingWindow(long index, long size) {
  return 0.54 - 0.46 * cos(2 * M_PI * index/double(size));
}

std::vector<double> SpectralProcessing::hammingWindow(std::vector<double> input) {
  std::vector<double> output(input.size());
  for (long i = 0; i < input.size(); i++)
    output[i] = input[i] * hammingWindow(i, input.size());
  return output;
}

double SpectralProcessing::quinnKappa(double in) {
  double firstTerm = log(3*in*in + 6*in + 1)/4.;

  double top = in + 1 - sqrt(2/3.);
  double bot = in + 1 + sqrt(2/3.);
  double secondTerm = sqrt(6)/24. * log(top/bot);

  return firstTerm - secondTerm;
}

double SpectralProcessing::quinnsFreqEstimator(std::complex<double> left,
                                               std::complex<double> mid,
                                               std::complex<double> right) {
  double betaM1 = (left/mid).real();
  double betaP1 = (right/mid).real();
  
  double deltaM1 = -betaM1/(betaM1 - 1);
  double deltaP1 = betaP1/(betaP1 - 1);

  double delta = (deltaM1 + deltaP1)/2.
                 + quinnKappa(deltaP1 * deltaP1)
                 - quinnKappa(deltaM1 * deltaM1);

  return delta;
}

std::complex<double> SpectralProcessing::quinnC(double delta, short indexK) {
  std::complex<double> I(0,1);
  std::complex<double> top = exp(2 * M_PI * I * delta) - 1.;
  std::complex<double> bot = 4 * M_PI * I * (delta - indexK);
  return top/bot;
}

double SpectralProcessing::quinnsAmpEstimator(double delta,
                          std::complex<double> left,
                          std::complex<double> mid,
                          std::complex<double> right) {

  std::complex<double> top = 0;
  top += left * std::conj(quinnC(delta, -1));
  top += mid * std::conj(quinnC(delta, 0));
  top += right * std::conj(quinnC(delta, 1));

  double bot = 0;
  for (long k = -1; k <= 1; k++) {
    double absC = std::abs(quinnC(delta, k));
    bot += absC * absC;
  }

  return std::abs(top)/bot;
}

std::complex<double> SpectralProcessing::fftToComplex(fftw_complex * fft, long bin) {
  return std::complex<double> (fft[bin][0], fft[bin][1]);
}

// In the future this should be normalized by the Equal loudness contour
// for accuracy to human hearing... 

std::vector< std::pair<double, double> > SpectralProcessing::findPeaks(
    fftw_complex * fft, long fftSize, long sampleRate) {
  std::vector< std::pair<double, double> > output;
  std::complex<double> left, mid, right;
  for (long bin = 1; bin < fftSize - 1; bin ++) {
    left = fftToComplex(fft, bin - 1);
    mid = fftToComplex(fft, bin);
    right = fftToComplex(fft, bin + 1);
    if (std::abs(mid) > std::abs(left) and std::abs(mid) > std::abs(right) ) {
      // A potential peak
      double delta = quinnsFreqEstimator(left, mid, right);
      if (fabs(delta) < 0.5) {
        // A definite peak
        double peakBin = bin + delta;
        double peakFreq = (peakBin * sampleRate)/(2*fftSize + 1);
        double peakAmp = quinnsAmpEstimator(delta, left, mid, right);
        std::pair<double, double> peak(peakFreq, peakAmp);
        output.push_back(peak);
      }
    }
  }
  return output;
}
