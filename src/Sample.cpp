#include <vector>
#include <iostream>

#include <SampleSorter/Sample.hpp>

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
  getWaves();
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

