#include <vector>

#include <fftw3.h>

#include "sample_sorter/beat.hpp"

using namespace sample_sorter;

Beat::Beat(
    const std::vector<std::vector<double>> & audio,
    long sample_rate) {
  // Compute the onset energy
  std::vector<std::vector<double>> onsets = onset_energy(audio);
}

Beat::onset_energy(const std::vector<std::vector<double>> & audio) {
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
  std::vector<double> secToLastMagnitude(fftSize);
  std::vector<double> lastPhase(fftSize);
  std::vector<double> secToLastPhase(fftSize);

  long maxHop = audio[0].size()/hopSize - windowRatio + 1;

  std::vector<double> output(maxHop);
  std::fill(output.begin(), output.end(), 0);

  for (long channel = 0; channel < audio.size(); channel ++) {
    std::fill(lastMagnitude.begin(), lastMagnitude.end(), 0);
    std::fill(secToLastMagnitude.begin(), secToLastMagnitude.end(), 0);
    std::fill(lastPhase.begin(), lastPhase.end(), 0);
    std::fill(secToLastPhase.begin(), secToLastPhase.end(), 0);

    for (long hop = 0; hop < maxHop; hop++) {
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
        std::complex<double> expected = (2*lastMagnitude[i] - secToLastMagnitude[i]) * std::exp(phaseShift);

        // Euclidean distance of expected and actual
        output[hop] += std::hypot(fft[i][0] - std::real(expected),
                                  fft[i][1] - std::imag(expected));
      }

      //store the current values
      for (long i = 0; i < fftSize; i++) {

}




#include <cmath>
#include <random>
#include <algorithm>


#include "SampleSorter/Tempo.hpp"
#include "SampleSorter/Units.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "Plotting/Plotting.hpp"

const long Tempo::HOP_SIZE = 1024;
const short Tempo::WINDOW_RATIO = 2;
// no tempos less than 1 beat per second
const double Tempo::HIGH_PASS = 1.;
// no tempos greater than 1000bpm
const double Tempo::LOW_PASS = 1000/60;

Tempo(
    const std::vector<std::vector<double> > & audio,
    long sampleRate) {

  std::vector<double> onsets =
    SpectralProcessing::onsetEnergy(audio, HOP_SIZE, windowRatio);
}


void Tempo::fineTuneTempo(
    const double percentageError,
    const int steps,
    const std::vector<double> & onsets
    ) {
  double guessTempo = tempo;

  std::vector<std::pair<double, double>> tempoValues(steps + 1);

  double minValue = onsets.size();
  double minTempo = guessTempo;
  int minOneBin = 0;

  // Guess a tempo in the range
  // iterate over bins, because
  // the error will be linear with number of bins
  double midBins = Units::tempoToBins(guessTempo, HOP_SIZE, sampleRate);
  double minBins = midBins * (1 - percentageError);
  double maxBins = midBins * (1 + percentageError);

  for (double bins = minBins;
      bins < maxBins;
      bins += (maxBins - minBins)/double(steps)) {

    tempo = Units::binsToTempo(bins, HOP_SIZE, sampleRate);

    // guess a one in the range
    for (theOneBin = 0;
        theOneBin < bins;
        theOneBin ++) {
        
      // calculate the value
      double newValue = getValue(onsets, true);

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
    const int steps) {
  // the first beat
  double beatOneBins = Units::tempoToBins(tempo, HOP_SIZE, sampleRate);

  double minValue = onsets.size();
  double minOneBin = 0;

  // for all possible ones
  for (theOneBin = 0;
      theOneBin < beatOneBins;
      theOneBin += beatOneBins/double(steps)) {
    double value = getValue(onsets, true);

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
    long sampleRate
    ) const {

  double distance = 
    (bin - theOneBin)/Units::tempoToBins(tempo, HOP_SIZE, sampleRate);
  if (bidirectional) {
    distance = std::abs(distance - std::round(distance));
  } else {
    distance = distance - std::floor(distance);
  }

  return distance;
}

double Tempo::getValue(
    const std::vector<double> & onsets,
    bool bidirectional
    ) const {

  double value = 0;
  for (size_t bin = 0; bin < onsets.size(); bin ++) {
    value += distanceFromBeat(bin, bidirectional, sampleRate) * onsets[bin];
  }
  
  return value;
}

Tempo::Tempo(
    const std::vector< std::vector<double> > & audio,
    long sampleRate_,
    double percentageError,
    int tempoSteps,
    int oneSteps) : sampleRate(sampleRate_) {

  long windowRatio = 2;
  std::vector<double> onsets =
    SpectralProcessing::onsetEnergy(audio, HOP_SIZE, windowRatio);

  findCorrelationTempo(onsets);

  fineTuneTempo(
      percentageError,
      tempoSteps,
      onsets
      );

  fineTuneTheOne(onsets, oneSteps);
}

void Tempo::findCorrelationTempo(
    const std::vector<double> & onsets
    ) {

  // filter out frequencies less than 2 per clip
  // Clips with tempo must contain at least 2 beats
  double highPass = 2./Units::binsToSeconds(onsets.size(), HOP_SIZE, sampleRate);
  highPass = (std::max)(highPass, HIGH_PASS);

  std::vector<double> correlation = 
    SpectralProcessing::autoCorrelation(
        onsets,
        sampleRate/HOP_SIZE,
        highPass,
        LOW_PASS);

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
    SpectralProcessing::findPeaks(fft, fftSize, sampleRate/HOP_SIZE);

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
    std::string title,
    std::string xaxis,
    std::string yaxis
    ) const {

  std::vector<double> beats;
  double beat = getTheOne();
  while (beat < Units::binsToSeconds(onsets.size(), HOP_SIZE, sampleRate)) {
    beats.push_back(beat);
    beat += Units::tempoToSeconds(tempo);
  }

  double binsToSeconds = Units::binsToSeconds(1, HOP_SIZE, sampleRate);
  Plotting::plotVector(onsets, title, xaxis, yaxis, binsToSeconds, 0, false, 0, beats);
}

