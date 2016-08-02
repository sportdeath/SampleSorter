#ifndef SPECTRAL_PROCESSING_H
#define SPECTRAL_PROCESSING_H

#include <complex>
#include <fftw3.h>

class SpectralProcessing {
  public:

    static double hammingWindow(long index, long size);

    static void hammingWindow(std::vector<double> & output, const std::vector<double> & input);

    static std::complex<double> fftToComplex(fftw_complex * fft, long bin);

    static double quinnKappa(double input);

    static double quinnsFreqEstimator(std::complex<double> left,
                                      std::complex<double> mid,
                                      std::complex<double> right);

    static std::complex<double> quinnC(double delta, short indexK);

    static double quinnsAmpEstimator(double delta,
                                     std::complex<double> left,
                                     std::complex<double> mid,
                                     std::complex<double> right);

    /**
     * Returns all peaks (frequency, value)
     * of a given Fourier transform
     */
    static std::vector< std::pair<double, double> > findPeaks(
        fftw_complex * fft, long fftSize, long sampleRate);

    static std::vector<double> onsetEnergy(
        std::vector< std::vector<double> > audio,
        long hopSize,
        long windowRatio);

    static std::vector<double> autoCorrelation(std::vector<double> signal,
                                               double sampleRate = 1,
                                               double highPass = 0,
                                               double lowPass = 1);

    static double getEnergy(const fftw_complex * fft, long fftSize);
    static double getAverageEnergyPerSample(const fftw_complex * fft, long fftSize);
    static double getAverageEnergyPerSecond(const fftw_complex * fft, long fftSize, long sampleRate);
};

#endif
