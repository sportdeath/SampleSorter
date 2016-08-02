#include <vector>
#include <cmath>
#include <iostream>

#include "SampleSorter/TimeDomainProcessing.hpp"

double TimeDomainProcessing::getEnergy(const std::vector<double> & signal) {
  double total = 0;
  for (long i = 0; i < signal.size(); i++) {
    total += signal[i] * signal[i];
  }
  return total;
}

double TimeDomainProcessing::getAverageEnergyPerSample(const std::vector<double> & signal) {
  return getEnergy(signal)/double(signal.size());
}

double TimeDomainProcessing::getAverageEnergyPerSecond(
    const std::vector<double> & signal, long sampleRate) {
  return getAverageEnergyPerSample(signal) * sampleRate;
}

void TimeDomainProcessing::normalizeEnergy(
    std::vector<std::vector<double> > & output,
    const std::vector<std::vector<double> > & audio,
    long sampleRate) {


  double averageEnergy = 0;
  for (long channel = 0; channel < audio.size(); channel ++) {
    averageEnergy += getAverageEnergyPerSecond(audio[channel], sampleRate);
  }
  averageEnergy /= double(audio.size());

  double rootEnergy = std::sqrt(averageEnergy);

  output.resize(audio.size());
  for (long channel = 0; channel < audio.size(); channel ++) {
    output[channel].resize(audio[channel].size());
    for (long i = 0; i < audio[channel].size(); i++) {
      output[channel][i] = audio[channel][i]/rootEnergy;
    }
  }
}
