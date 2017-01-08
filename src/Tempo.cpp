#include <vector>
#include <cmath>
#include <random>

#include <fftw3.h>

#include "SampleSorter/Tempo.hpp"
#include "SampleSorter/Units.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "Plotting/Plotting.hpp"

Tempo::Tempo() {
  tempo = 1;
  theOneBin = 0;
}

double Tempo::getTheOneBin() const {
  return theOneBin;
}

double Tempo::getTheOne(
    long hopSize, 
    long sampleRate
    ) const {
  return Units::binsToSeconds(theOneBin, hopSize, sampleRate);
}

Tempo::Tempo(double tempo_, double theOneBin_) {
  tempo = tempo_;
  theOneBin = theOneBin_;
}

void Tempo::fineTuneTempo(
    const double percentageError,
    const int steps,
    const std::vector<double> & onsets,
    const long hopSize,
    const long sampleRate
    ) {
  double guessTempo = tempo;

  double minValue = onsets.size();
  double minTempo = guessTempo;
  int minOneBin = 0;

  // Guess a tempo in the range
  // iterate over bins, because
  // the error will be linear with number of bins
  double midBins = Units::tempoToBins(guessTempo, hopSize, sampleRate);
  double minBins = midBins * (1 - percentageError);
  double maxBins = midBins * (1 + percentageError);
  for (double bins = minBins;
      bins < maxBins;
      bins += (maxBins - minBins)/double(steps)) {

    tempo = Units::binsToTempo(bins, hopSize, sampleRate);
    
    // guess a one in the range
    for (theOneBin = 0;
        theOneBin < bins;
        theOneBin ++) {
      // calculate the value
      double newValue = getValue(onsets, true, hopSize, sampleRate);

      // choose one that minimize energy off of beat
      if (newValue < minValue) {
        minValue = newValue;
        minTempo = tempo;
        minOneBin = theOneBin;
      }
    }

  }

  tempo = minTempo;
  theOneBin = minOneBin;
}

void Tempo::fineTuneTheOne(
    const std::vector<double> & onsets,
    const int steps,
    const long hopSize,
    const long sampleRate) {
  // the first beat
  double beatOneBins = Units::tempoToBins(tempo, hopSize, sampleRate);

  double minValue = onsets.size();
  double minOneBin = 0;

  // for all possible ones
  for (theOneBin = 0;
      theOneBin < beatOneBins;
      theOneBin += beatOneBins/double(steps)) {
    double value = getValue(onsets, false, hopSize, sampleRate);

    if (value < minValue) {
      minValue = value;
      minOneBin = theOneBin;
    }
  }
  theOneBin = minOneBin;
}

double Tempo::getTempo() const {
  return tempo;
}

double Tempo::distanceFromBeat(
    double bin,
    bool bidirectional,
    long hopSize,
    long sampleRate
    ) const {
  double distance = 
    (bin - theOneBin)/Units::tempoToBins(tempo, hopSize, sampleRate);
  if (bidirectional) {
    distance = std::abs(distance - std::round(distance));
  } else {
    distance = distance - std::floor(distance);
  }
  return distance;
}

double Tempo::getValue(
    const std::vector<double> & onsets,
    bool bidirectional,
    long hopSize,
    long sampleRate
    ) const {

  double value = 0;
  for (size_t bin = 0; bin < onsets.size(); bin ++) {
    value += distanceFromBeat(bin, hopSize, sampleRate, bidirectional) * onsets[bin];
  }
  
  return value;
}

Tempo::Tempo(
    const std::vector< std::vector<double> > & audio,
    long sampleRate,
    double percentageError,
    int tempoSteps,
    int oneSteps) {

  std::cout << audio.size() << ", " << audio[0].size() << std::endl;
  std::cout << sampleRate << std::endl;

  //std::cout << "starting tempo" << std::endl;

  long hopSize = 1024;
  long windowRatio = 2;
  std::vector<double> onsets =
    SpectralProcessing::onsetEnergy(audio, hopSize, windowRatio);

  findCorrelationTempo(onsets, hopSize, sampleRate);

  std::cout << "guessed tempo: " << getTempo() * 60. << std::endl;

  fineTuneTempo(
      percentageError,
      tempoSteps,
      onsets,
      hopSize,
      sampleRate
      );

  std::cout << "fine tuned tempo: " << getTempo() * 60. << std::endl;

  //aCoefficients.resize(degrees);
  //bCoefficients.resize(degrees);
  //std::fill(aCoefficients.begin(), aCoefficients.end(), 0);
  //std::fill(bCoefficients.begin(), bCoefficients.end(), 0);
  //gradientDescent(onsets, hopSize, sampleRate);

  std::cout << "the One: " << getTheOne(hopSize, sampleRate) << std::endl;

  fineTuneTheOne(onsets, oneSteps, hopSize, sampleRate);

  std::cout << "the fine tuned One: " << getTheOne(hopSize, sampleRate) << std::endl;

  //for (long i = 0; i < getDegree(); i++) {
    //std::cout << "a[" << i << "]=" << aCoefficients[i] << std::endl;
    //std::cout << "a[" << i << "]=" << bCoefficients[i] << std::endl;
  //}

  plotOnsetsWithBeats(onsets, hopSize, sampleRate);
  //plotBeats(0.001);
  //plotTempo(0.001);
  //

}

void Tempo::findCorrelationTempo(
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

  std::vector<double> correlation = 
    SpectralProcessing::autoCorrelation(
        onsets,
        sampleRate/hopSize,
        highPass,
        lowPass);

  // only take first half of correlation
  long fftSize = correlation.size()/4 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(
      correlation.size()/2,
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

  Plotting::plotPair(peaks);

  fftw_free(fft);

  double max = 0;
  tempo = 0;
  for (long i = 0; i < peaks.size(); i++) {
    if (peaks[i].second > max) {
      max = peaks[i].second;
      tempo = peaks[i].first;
    }
  }
}

// plot onsets with beats
void Tempo::plotOnsetsWithBeats(
    const std::vector<double> & onsets,
    long hopSize,
    long sampleRate
    ) const {

  std::vector<double> beats;
  double beat = getTheOne(hopSize, sampleRate);
  while (beat < Units::binsToSeconds(onsets.size(), hopSize, sampleRate)) {
    beats.push_back(beat);
    beat += Units::tempoToSeconds(tempo);
  }

  std::vector<std::pair<double, double> > onsetSeconds(onsets.size());

  for (long bin = 0; bin < onsets.size(); bin++) {
    onsetSeconds[bin] = std::make_pair(Units::binsToSeconds(bin, hopSize, sampleRate), onsets[bin]);
  }

  Plotting::plotLineAndMarkers(onsetSeconds, beats, 0.5);
}

