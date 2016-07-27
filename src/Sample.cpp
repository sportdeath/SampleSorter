#include <vector>
#include <iostream>

#include <SampleSorter/Sample.hpp>
#include <SampleSorter/Octave.hpp>
#include <SampleSorter/SpectralProcessing.hpp>

#include "gnuplot-iostream.h"

Sample::Sample(std::string file_) {
  file = file_;
}

void Sample::process() {
  tune();
  findBeat();
  //findChords();
}


std::string Sample::getFile() {
  return file;
}

void Sample::tune() {
  std::vector<std::vector<double> > a = getWaves();
  Octave o(a, 1200, getSampleRate());
  long cents = o.tune();
  std::cout << "tuned by " << cents << " cents" << std::endl;
  return;
}

void Sample::findBeat() {
  std::vector<std::vector<double> > a = getWaves();
  long hopsize = 128;

  double tempo 
    = SpectralProcessing::tempoDetection(a, hopsize, 2, getSampleRate());

  std::cout << "tempo: " << tempo*60 << " bpm" << std::endl;
  return;
}

//void Sample::findChords() {
  //return;
//}

bool Sample::isHarmonic() {
  return isHarmonic_;
}

bool Sample::hasBeat() {
  return hasBeat_;
}

double Sample::getBeat() {
  return beat;
}

//std::vector<Chord> Sample::getChords() {
  //return chords;
//}
