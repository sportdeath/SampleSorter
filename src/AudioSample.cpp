#include <vector>
#include <iostream>

#include <fftw3.h>

#include "SampleSorter/AudioSample.hpp"
#include "SampleSorter/TimeDomainProcessing.hpp"
#include "SampleSorter/SpectralProcessing.hpp"
#include "SampleSorter/EqualLoudness.hpp"
#include "SampleSorter/Octave.hpp"
#include "SampleSorter/Units.hpp"
#include "SampleSorter/Tempo.hpp"

#include "Plotting/Plotting.hpp"

AudioSample::AudioSample(
    std::vector<std::vector<double> > & audio,
    long sampleRate) {

  totalSeconds = Units::samplesToSeconds(audio[0].size(), sampleRate);

  // Filter the audio to get the sound
  // as heard by the human ear
  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);

  // Normalize the audio to have unit energy per second
  // TODO

  // Tune the audio
  tune(filteredAudio, sampleRate);
  // Find the tempo
  findTempo(filteredAudio, sampleRate);
  // Find out what chords dominate the audio
  findChords(filteredAudio, sampleRate);
}

void AudioSample::tune(
    std::vector<std::vector<double> > & audio,
    long sampleRate) {

  // Put the audio into an octave format
  Octave oct(audio, 1200., sampleRate);

  // Tune it to the nearest cent
  tuningCents = oct.tune();
}

void AudioSample::findTempo(
    std::vector<std::vector<double> > & audio,
    long sampleRate) {

  // Compute the tempo
  tempo = Tempo(audio, sampleRate);
}

void AudioSample::findChords(
    std::vector<std::vector<double> > & audio,
    long sampleRate) {

  // Get how many samples are in each window
  size_t windowSize = Units::tempoToSamples(tempo.getTempo(), sampleRate);
  // Determine where the start of the sample is
  size_t startSample = Units::secondsToSamples(tempo.getTheOne(), sampleRate);

  int numChords = std::floor((totalSeconds - tempo.getTheOne()) * tempo.getTempo());
  chords.resize(numChords);

  std::vector<std::vector<double>> window(audio.size(), std::vector<double>(windowSize));

  for (auto it = chords.begin(); it != chords.end(); it++) {

    // Copy the content into the window
    for (size_t channel = 0; channel < audio.size(); channel ++) {
      for (size_t i = 0; i < windowSize; i++) {
        // Apply a hamming window
        window[channel][i] = 
          audio[channel][i + startSample] *
          SpectralProcessing::hammingWindow(i, windowSize);
      }
    }

    // Add the octave
    chords.emplace(window, 12, sampleRate, tuningCents);
    
    // Move to the next window
    startSample += windowSize;
  }
}
