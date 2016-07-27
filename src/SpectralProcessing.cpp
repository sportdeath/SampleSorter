#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

#include "SampleSorter/SpectralProcessing.hpp"

double SpectralProcessing::hammingWindow(long index, long size) {
  return 0.54 - 0.46 * cos(2 * M_PI * index/double(size));
}

void SpectralProcessing::hammingWindow(std::vector<double> & output, const std::vector<double> & input) {
  for (long i = 0; i < input.size(); i++)
    output[i] = input[i] * hammingWindow(i, input.size());
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


// From Bello et al.
// On the use of phase and energy for musical onset detection in the complex domain
std::vector<double> SpectralProcessing::onsetEnergy(
    std::vector< std::vector<double> > audio, 
    long hopSize, 
    long windowRatio) {

  long windowSize = hopSize * windowRatio;

  long fftSize = windowSize/2 + 1;
  std::vector<double> window(windowSize);
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
                                              window.data(),
                                              fft,
                                              FFTW_MEASURE);

  std::vector<double> lastMagnitude(fftSize);
  std::vector<double> lastPhase(fftSize);
  std::vector<double> secToLastPhase(fftSize);

  std::vector<double> output(audio[0].size()/hopSize);
  std::fill(output.begin(), output.end(), 0);

  for (long channel = 0; channel < audio.size(); channel ++) {
    std::fill(lastMagnitude.begin(), lastMagnitude.end(), 0);
    std::fill(lastPhase.begin(), lastPhase.end(), 0);
    std::fill(secToLastPhase.begin(), secToLastPhase.end(), 0);

    for (long hop = 0; hop < audio[channel].size()/hopSize; hop++) {
      // fill window with windowSize values
      // beginning at hopSize * hop
      // using hamming window
      for (long i = 0; i < windowSize; i ++) {
        window[i] = hammingWindow(i, windowSize)
                    * audio[channel][hop * hopSize + i];
      }

      fftw_execute(fftPlan);

      for (long i = 0; i < fftSize; i++) {
        //calculate expected values
        std::complex<double> phaseShift(0, 2*lastPhase[i] - secToLastPhase[i]);
        std::complex<double> expected = lastMagnitude[i] * std::exp(phaseShift);

        // Euclidean distance of expected and actual
        output[hop] += std::hypot(fft[i][0] - std::real(expected),
                                  fft[i][1] - std::imag(expected));
      }

      //store the current values
      for (long i = 0; i < fftSize; i++) {
        lastMagnitude[i] = std::hypot(fft[i][0], fft[i][1]);
        secToLastPhase[i] = lastPhase[i];
        lastPhase[i] = std::atan2(fft[i][0], fft[i][1]);
      }
    }
  }

  for (long i = 0; i < output.size(); i++) {
    output[i] = output[i]/windowSize;
  }

  return output;
}

std::vector<double> SpectralProcessing::autoCorrelation(std::vector<double> signal) {
  // Should this be windowed??

  // FFT
  long fftSize = signal.size()/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlanFor = fftw_plan_dft_r2c_1d(signal.size(),
                                              signal.data(),
                                              fft,
                                              FFTW_ESTIMATE);

  fftw_execute(fftPlanFor);
  fftw_destroy_plan(fftPlanFor);

  // multiply with complex conjugate -> (a+bi)(a-bi) = (a^2 + b^2 + 0i)
  for (long i = 0; i < fftSize; i++) {
    fft[i][0] = fft[i][0] * fft[i][0] + fft[i][1] * fft[i][1];
    fft[i][1] = 0;
  }

  std::vector<double> output(signal.size());

  // Inverse FFT
  fftw_plan fftPlanInv = fftw_plan_dft_c2r_1d(output.size(),
                                              fft,
                                              output.data(),
                                              FFTW_ESTIMATE);

  fftw_execute(fftPlanInv);
  fftw_destroy_plan(fftPlanInv);
  fftw_free(fft);

  return output;
}


double SpectralProcessing::tempoDetection(
    std::vector< std::vector<double> > audio, 
    long hopSize, 
    long windowRatio,
    long sampleRate) {


  std::vector<double> onsets = onsetEnergy(audio, hopSize, windowRatio);
  std::vector<double> correlation = autoCorrelation(onsets);

  long fftSize = correlation.size()/2 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(correlation.size(),
                                              correlation.data(),
                                              fft,
                                              FFTW_ESTIMATE);

  fftw_execute(fftPlan);
  fftw_destroy_plan(fftPlan);

  std::vector<std::pair<double, double> > peaks = 
    findPeaks(fft, fftSize, sampleRate/hopSize);

  fftw_free(fft);

  double max = 0;
  double maxFreq = 0;
  for (long i = 0; i < peaks.size(); i++) {
    if (peaks[i].second > max) {
        //and peaks[i].first * 60 > 20) {
      max = peaks[i].second;
      maxFreq = peaks[i].first;
    }
  }

  return maxFreq;
}


// get onset energy
// get autocorrelation of it
// then take the first half of it and compute the fourier transform
// 
// onset energy has samplingRate/hopSize samples per second
// autocorrelation has same sampling rate
// transform of autocorrelation...
