#include <vector>
#include <cmath>
#include <random>

#include <fftw3.h>

#include "SampleSorter/TempoFunction.hpp"
#include "SampleSorter/Units.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "Plotting/Plotting.hpp"

TempoFunction::TempoFunction() {
  avgTempo = 1;
  totalSeconds = 0;
}

double TempoFunction::getTotalSeconds() const {
  return totalSeconds;
}

double TempoFunction::getTheOne() const {
  return aCoefficients[0]/2.;
}

TempoFunction::TempoFunction(double avgTempo_,
              double theOne,
              double totalSeconds_) {
  aCoefficients.resize(1);
  bCoefficients.resize(1);
  avgTempo = avgTempo_;
  aCoefficients[0] = theOne * 2.;
  totalSeconds = totalSeconds_;
}

void TempoFunction::fineTuneAvgTempo(
    double guessTempo,
    double percentageError,
    double percentageStep,
    const std::vector<double> & onsets,
    long hopSize,
    long sampleRate
    ) {

  double minValue = onsets.size();
  double minTempo = guessTempo;
  double minOne = 0;

  for (double trialTempo = guessTempo/2.;
       trialTempo < guessTempo * 2.;
       trialTempo *= std::pow(2, 0.001)) {

    avgTempo = trialTempo;
    aCoefficients[0] = 0;

    double theFirstBeat = getKthBeat(1);

    for (double theOne = 0;
        theOne < Units::secondsToBins(theFirstBeat, hopSize, sampleRate);
        theOne ++) {

      aCoefficients[0] = theOne;

      double newValue = getValue(onsets, hopSize, sampleRate);
      if (newValue < minValue) {
        minValue = newValue;
        minTempo = trialTempo;
        minOne = theOne;
      }
    }
  }

  avgTempo = minTempo;
  aCoefficients[0] = minOne;
}

void TempoFunction::fineTuneTheOne(
    const std::vector<double> & onsets,
    long hopSize,
    long sampleRate) {
  // set the one to zero
  aCoefficients[0] = 0;

  // get the first beat
  double beatOne = getKthBeat(1);

  double minValue = onsets.size(), minOne = 0;

  // for all possible ones
  for (double onePosition = 0;
      onePosition < beatOne;
      onePosition += 0.001) {
    aCoefficients[0] = onePosition;
    double value = 0;
    for (long bin = 0; bin < onsets.size(); bin ++) {
      double binSeconds = Units::binsToSeconds(bin, hopSize, sampleRate);
      value += onsets[bin] * std::fmod(getBeatNum(binSeconds), 1.);
    }
    if (value < minValue) {
      minValue = value;
      minOne = onePosition;
    }
  }
  aCoefficients[0] = minOne;
}

double TempoFunction::getAvgTempo() const {
  return avgTempo;
}

size_t TempoFunction::getDegree() const {
  return aCoefficients.size();
}

double TempoFunction::getBeatNum(double time) const {
  // actual tempo
  double beatNum = avgTempo * time;
  // constant offset
  beatNum += aCoefficients[0]/2.;

  // Fourier series approximation
  for (size_t i = 1; i < getDegree(); i++) {
    beatNum += aCoefficients[i] * 
               std::cos((2.*M_PI*i*time)/totalSeconds);
    beatNum += bCoefficients[i] * 
               std::sin((2.*M_PI*i*time)/totalSeconds);
  }

  return beatNum;
}

// The derivative
double TempoFunction::getTempo(double time) const {
  double tempo = avgTempo;

  // Fourier series approximation
  for (size_t i = 1; i < getDegree(); i++) {
    tempo -= aCoefficients[i] * 
             ((2*M_PI*i)/totalSeconds) *
             std::sin(2.*M_PI*i*time/totalSeconds);

    tempo += bCoefficients[i] * 
             ((2*M_PI*i)/totalSeconds) *
             std::cos(2.*M_PI*i*time/totalSeconds);
  }

  return tempo;
}

double TempoFunction::distanceFromBeat(double time) const {
  return 1 - std::cos(2.*M_PI*getBeatNum(time));
}

double TempoFunction::getValue(
    const std::vector<double> & onsets,
    long hopSize,
    long sampleRate
    ) const {
  double value = 0;
  for (size_t bin = 0; bin < onsets.size(); bin ++) {
    double time = Units::binsToSeconds(bin, hopSize, sampleRate);
    value += distanceFromBeat(time) * onsets[bin];
  }
  
  return value;
}

void TempoFunction::getBinGradient(
    size_t bin,
    std::vector<double> & aNabla,
    std::vector<double> & bNabla,
    const std::vector<double> & onsets,
    long hopSize,
    long sampleRate
    ) const {

  double time = Units::binsToSeconds(bin,hopSize, sampleRate);

  double multiplier = onsets[bin];
  multiplier *= 2.*M_PI * std::sin(2.*M_PI*getBeatNum(time));

  aNabla[0] = 0.5 * multiplier;

  for (size_t i = 1; i < getDegree(); i++) {
    aNabla[i] = std::cos(2.*M_PI*i*time/totalSeconds);
    aNabla[i] *= multiplier;

    bNabla[i] = std::sin(2.*M_PI*i*time/totalSeconds);
    bNabla[i] *= multiplier;
  }
}

void TempoFunction::gradientDescent(
    const std::vector<double> & onsets,
    long hopSize, 
    long sampleRate
    ) {

  std::vector<double> aNabla(getDegree());
  std::vector<double> bNabla(getDegree());

  std::vector<double> aDelta(getDegree()); 
  std::fill(aDelta.begin(), aDelta.end(), 0);
  std::vector<double> bDelta(getDegree()); 
  std::fill(bDelta.begin(), bDelta.end(), 0);

  double learningRate;
  double mass = 0.5;

  std::srand(std::time(0));

  long numIterations = onsets.size() * 4;
  //std::cout << "numIterations" << numIterations << std::endl;
  std::vector<double> values(numIterations);

  for (long trial = 0; trial < numIterations; trial++) {
    learningRate = 0.000001*(numIterations - trial)/double(numIterations);

    // choose a random bin
    size_t bin = std::rand() % onsets.size();

    getBinGradient(bin, aNabla, bNabla, onsets, hopSize, sampleRate);

    for (size_t i = 0; i < getDegree(); i++) {
      aDelta[i] = learningRate * aNabla[i] + mass * aDelta[i];
      aCoefficients[i] = aCoefficients[i] - aDelta[i];

      bDelta[i] = learningRate * bNabla[i] + mass * bDelta[i];
      bCoefficients[i] = bCoefficients[i] - bDelta[i];
    }

    values[trial] = getValue(onsets, hopSize, sampleRate);
  }

}


double TempoFunction::getKthBeat(double k, double accuracy) const {
  //std::cout << "starting" << k << std::endl;
  double guess = k/avgTempo + aCoefficients[0];

  double error = getBeatNum(guess) - k;
  while (std::abs(error) > accuracy) {
    guess -= error/getTempo(guess);
    //std::cout << "tempo: " << getTempo(guess) << std::endl;
    //std::cout << "error: " << error << std::endl;
    //std::cout << "guess: " << guess << std::endl;
    error = getBeatNum(guess) - k;
  }

  //std::cout << "ending" << k << std::endl;

  return guess;
}

TempoFunction::TempoFunction(
    std::vector< std::vector<double> > & audio,
    long sampleRate,
    size_t degrees, 
    double percentageError,
    double percentageStep) {

  std::cout << audio.size() << ", " << audio[0].size() << std::endl;
  std::cout << sampleRate << std::endl;

  //std::cout << "starting tempo" << std::endl;

  totalSeconds = Units::samplesToSeconds(audio[0].size(), sampleRate);

  long hopSize = 1024;
  long windowRatio = 2;
  std::vector<double> onsets =
    SpectralProcessing::onsetEnergy(audio, hopSize, windowRatio);

  double guessTempo = correlationTempo(onsets, hopSize, sampleRate);

  std::cout << "guessed tempo: " << guessTempo * 60. << std::endl;

  aCoefficients.resize(1);
  bCoefficients.resize(1);
  std::fill(aCoefficients.begin(), aCoefficients.end(), 0);
  std::fill(bCoefficients.begin(), bCoefficients.end(), 0);

  fineTuneAvgTempo(guessTempo, 
                   percentageError, 
                   percentageStep, 
                   onsets, 
                   hopSize, 
                   sampleRate);

  std::cout << "fine tuned tempo: " << getAvgTempo() * 60. << std::endl;

  //aCoefficients.resize(degrees);
  //bCoefficients.resize(degrees);
  //std::fill(aCoefficients.begin(), aCoefficients.end(), 0);
  //std::fill(bCoefficients.begin(), bCoefficients.end(), 0);
  //gradientDescent(onsets, hopSize, sampleRate);

  fineTuneTheOne(onsets, hopSize, sampleRate);

  //for (long i = 0; i < getDegree(); i++) {
    //std::cout << "a[" << i << "]=" << aCoefficients[i] << std::endl;
    //std::cout << "a[" << i << "]=" << bCoefficients[i] << std::endl;
  //}

  //plotOnsetsWithBeats(onsets, hopSize, sampleRate);
  //plotBeats(0.001);
  //plotTempo(0.001);
  //

}

double TempoFunction::correlationTempo(
                     const std::vector<double> & onsets, 
                     long hopSize, 
                     long sampleRate) {

  // filter out frequencies less than 2 per clip
  // Clips with tempo must contain at least 2 beats
  double highPass = 2./Units::binsToSeconds(onsets.size(), hopSize, sampleRate);
  // Also beat cannot be below 1/8 persecond
  highPass = std::max(highPass, 1.);

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
    correlation[i] = correlation[i] 
              * SpectralProcessing::hammingWindow(i, correlation.size()/2);
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

  return tempo;
}

void TempoFunction::plotTempo(double resolution) const {
  std::vector<double> tempo(totalSeconds/resolution);
  for (long i = 0; i*resolution < totalSeconds; i++) {
    tempo[i] = 60.*getTempo(i*resolution);
  }
  Plotting::plotVector(tempo);
}

void TempoFunction::plotBeats(double resolution) const {
  std::vector<double> beats(totalSeconds/resolution);
  for (long i = 0; i*resolution < totalSeconds; i++) {
    beats[i] = getBeatNum(i*resolution);
  }
  Plotting::plotVector(beats);
}

void TempoFunction::plotOnsetsWithBeats(
    const std::vector<double> & onsets,
    long hopSize,
    long sampleRate
    ) const {

  std::vector<double> beats;
  long k = 0;
  double beat;
  do {
    beat = getKthBeat(k);
    beats.push_back(beat);
    k += 1;
  } while (beat < totalSeconds);

  std::vector<std::pair<double, double> > onsetSeconds(onsets.size());

  for (long bin = 0; bin < onsets.size(); bin++) {
    onsetSeconds[bin] = std::make_pair(Units::binsToSeconds(bin, hopSize, sampleRate), onsets[bin]);
  }

  Plotting::plotLineAndMarkers(onsetSeconds, beats, 0.5);
}

// plot onsets with beats
