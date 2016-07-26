#include <complex>
#include <fftw3.h>

class SpectralProcessing {
  public:

    static double hammingWindow(long index, long size);

    static std::vector<double> hammingWindow(std::vector<double> input);

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
};
