#include <vector>
#include <iostream>

#include <SampleSorter/Sample.hpp>
#include <SampleSorter/Octave.hpp>

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
  Octave o(a, 60, getSampleRate());
  o.plot();
  o.tune();
  // take fourier transform
  // wrap it to an octave
  // want to rotate it so the spectral energy
  // is as close to bins as possible.
  // minimize sum [ value@pos * distToNearestBin(pos) ]


  return;
}

void Sample::findBeat() {
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
