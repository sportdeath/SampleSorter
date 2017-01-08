#include <vector>
#include <iostream>

#include "SampleSorter/Units.hpp"

double Units::tempoToSeconds(double tempo) {
  return 1/tempo;
}

double Units::tempoToSamples(double tempo, long sampleRate) {
  return secondsToSamples(tempoToSeconds(tempo), sampleRate);
}

double Units::secondsToSamples(double seconds, long sampleRate) {
  return seconds * sampleRate;
}

double Units::secondsToTempo(double seconds) {
  return 1/seconds;
}

double Units::samplesToSeconds(double samples, long sampleRate) {
  return samples/double(sampleRate);
}

double Units::tempoToBins(double tempo, long hopSize, long sampleRate) {
  return sampleRate/(tempo * hopSize);
}

double Units::binsToSamples(double bin, long hopSize) {
  return bin * hopSize;
}

double Units::binsToSeconds(double bin, long hopSize, long sampleRate) {
  return samplesToSeconds(binsToSamples(bin, hopSize), sampleRate);
}

double Units::binsToTempo(double bin, long hopSize, long sampleRate) {
  return secondsToTempo(binsToSeconds(bin, hopSize, sampleRate));
}

double Units::samplesToBins(double samples, long hopSize) {
  return samples/double(hopSize);
}

double Units::secondsToBins(double seconds, long hopSize, long sampleRate) {
  return samplesToBins(secondsToSamples(seconds, sampleRate), hopSize);
}
