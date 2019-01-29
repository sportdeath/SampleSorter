#include <vector>

#include "sample_sorter/sample.hpp"
#include "sample_sorter/units.hpp"
#include "sample_sorter/octave.hpp"
#include "sample_sorter/beat.hpp"

using namespace sample_sorter;

Sample::Sample(
    std::vector<std::vector<double> > & audio,
    double sample_rate) {

  // Determine how long the sample is
  total_seconds_ = Units::samples_to_seconds(audio[0].size(), sample_rate);

  // Filter the audio to get the sound
  // as heard by the human ear
  // TODO Add in equal loudness
  //std::vector<std::vector<double>> filtered_audio =
    //EqualLoudness::filter(audio, sample_rate, filtered_audio);

  // Tune the audio
  tune(audio, sample_rate);
  // Find the tempo
  find_tempo(audio, sample_rate);
  // Find out what chords dominate the audio
  find_chords(audio, sample_rate);
}

void Sample::tune(
    std::vector<std::vector<double> > & audio,
    double sample_rate) {

  // Put the audio into an octave format
  Octave octave(
      audio, 
      Octave::CENTS_PER_OCTAVE,
      sample_rate);

  // Tune it to the nearest cent
  tuning_cents_ = octave.tune();
}

void Sample::find_tempo(
    std::vector<std::vector<double> > & audio,
    double sample_rate) {

  // Compute the tempo and the one
  beat = Beat(audio, sample_rate);
}

void Sample::find_chords(
    std::vector<std::vector<double> > & audio,
    double sample_rate) {

  // Get how many samples are in each window
  size_t window_size = Units::tempo_to_samples(beat.tempo(), sample_rate);
  std::vector<std::vector<double>> window(
      audio.size(), std::vector<double>(window_size));

  // Determine the number of chords there will be
  int num_chords = std::floor((total_seconds() - beat.the_one()) * beat.tempo());
  chords_.resize(num_chords);

  // Determine where the start of the chords is
  size_t start_sample = Units::seconds_to_samples(beat.the_one(), sample_rate);

  for (auto it = chords_.begin(); it != chords_.end(); it++) {

    // Copy the content into the window
    for (size_t channel = 0; channel < audio.size(); channel ++) {
      for (size_t i = 0; i < window_size; i++) {
        // Apply a hamming window
        window[channel][i] = 
          audio[channel][i + start_sample];
          // TODO Add in windowing
          //*
          //SpectralProcessing::hammingWindow(i, windowSize);
      }
    }

    // Compute the octave
    chords_.emplace(
        it, window, 
        Octave::SEMITONES_PER_OCTAVE, 
        sample_rate, 
        tuning_cents());
    
    // Move to the next window
    start_sample += window_size;
  }
}
