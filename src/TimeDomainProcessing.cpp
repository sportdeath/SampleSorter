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

double TimeDomainProcessing::getAverageEnergyPerBeat(
    const std::vector<double> & signal, double tempo, long sampleRate) {
  return getAverageEnergyPerSecond(signal, sampleRate)/tempo;
}

void TimeDomainProcessing::normalizeByEnergy(
    std::vector<std::vector<double> > & output,
    const std::vector<std::vector<double> > & audio,
    double energy) {
  double rootEnergy = std::sqrt(energy);

  output.resize(audio.size());
  for (long channel = 0; channel < audio.size(); channel ++) {
    output[channel].resize(audio[channel].size());
    for (long i = 0; i < audio[channel].size(); i++) {
      output[channel][i] = audio[channel][i]/rootEnergy;
    }
  }
}

void TimeDomainProcessing::unitEnergyPerSecond(
    std::vector<std::vector<double> > & output,
    const std::vector<std::vector<double> > & audio,
    long sampleRate) {
  double averageEnergy = 0;
  for (long channel = 0; channel < audio.size(); channel ++) {
    averageEnergy += getAverageEnergyPerSecond(audio[channel], sampleRate);
  }
  averageEnergy /= double(audio.size());
  normalizeByEnergy(output, audio, averageEnergy);
}

void TimeDomainProcessing::unitEnergyPerBeat(
    std::vector<std::vector<double> > & output,
    const std::vector<std::vector<double> > & audio,
    double tempo,
    long sampleRate) {
  double averageEnergy = 0;
  for (long channel = 0; channel < audio.size(); channel ++) {
    averageEnergy += getAverageEnergyPerBeat(audio[channel],tempo, sampleRate);
  }
  averageEnergy /= double(audio.size());
  normalizeByEnergy(output, audio, averageEnergy);
}
