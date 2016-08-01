#include <vector>
#include <string>

#include "SampleSorter/Octave.hpp"
#include "SampleSorter/MusicTheory.hpp"

const std::vector<std::string> MusicTheory::noteNames = 
    {"A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"};

const Octave MusicTheory::unison       (std::vector<double>({1,0,0,0,0,0,0,0,0,0,0,0}));
const Octave MusicTheory::minorSecond (std::vector<double>({1,1,0,0,0,0,0,0,0,0,0,0}));
const Octave MusicTheory::majorSecond (std::vector<double>({1,0,1,0,0,0,0,0,0,0,0,0}));
const Octave MusicTheory::minorThird  (std::vector<double>({1,0,0,1,0,0,0,0,0,0,0,0}));
const Octave MusicTheory::majorThird  (std::vector<double>({1,0,0,0,1,0,0,0,0,0,0,0}));
const Octave MusicTheory::perfectFifth(std::vector<double>({1,0,0,0,0,0,0,1,0,0,0,0}));

const Octave MusicTheory::majorTriad     (std::vector<double>({1,0,0,0,1,0,0,1,0,0,0,0}));
const Octave MusicTheory::minorTriad     (std::vector<double>({1,0,0,1,0,0,0,1,0,0,0,0}));
const Octave MusicTheory::diminishedTriad(std::vector<double>({1,0,0,1,0,0,1,0,0,0,0,0}));
const Octave MusicTheory::augmentedTriad (std::vector<double>({1,0,0,0,1,0,0,0,1,0,0,0}));

const Octave MusicTheory::sus2(std::vector<double>({1,0,1,0,0,0,0,1,0,0,0,0}));
const Octave MusicTheory::sus4(std::vector<double>({1,0,0,0,0,1,0,1,0,0,0,0}));

const Octave MusicTheory::minorSeventh(std::vector<double>({1,0,0,1,0,0,0,1,0,0,1,0}));
const Octave MusicTheory::minorMajorSeventh(std::vector<double>({1,0,0,1,0,0,0,1,0,0,0,1}));
const Octave MusicTheory::dominantSeventh(std::vector<double>({1,0,0,0,1,0,0,1,0,0,1,0}));
const Octave MusicTheory::majorSeventh(std::vector<double>({1,0,0,0,1,0,0,1,0,0,0,1}));
const Octave MusicTheory::fullDiminishedSeventh(std::vector<double>({1,0,0,1,0,0,1,0,0,1,0,0}));
const Octave MusicTheory::halfDiminishedSeventh(std::vector<double>({1,0,0,1,0,0,1,0,0,0,1,0}));
const Octave MusicTheory::augmentedSeventh(std::vector<double>({1,0,0,0,1,0,0,0,1,0,1,0}));
const Octave MusicTheory::augmentedMajorSeventh(std::vector<double>({1,0,0,0,1,0,0,0,1,0,0,1}));


const Octave MusicTheory::thirdlessDominantSeventh(std::vector<double>({1,0,0,0,0,0,0,1,0,0,1,0}));
const Octave MusicTheory::thirdlessMajorSeventh(std::vector<double>({1,0,0,0,0,0,0,1,0,0,0,1}));

const Octave MusicTheory::fifthlessMinorSeventh(std::vector<double>({1,0,0,1,0,0,0,0,0,0,1,0}));
const Octave MusicTheory::fifthlessMinorMajorSeventh(std::vector<double>({1,0,0,1,0,0,0,1,0,0,0,1}));
const Octave MusicTheory::fifthlessDominantSeventh(std::vector<double>({1,0,0,0,1,0,0,1,0,0,1,0}));
const Octave MusicTheory::fifthlessMajorSeventh(std::vector<double>({1,0,0,0,1,0,0,0,0,0,0,1}));

const std::vector<std::pair<std::string, Octave> > MusicTheory::chords = [] {
  std::vector<std::pair<std::string, Octave> > c;
  // intervals
  c.push_back(std::make_pair("Unison", unison));
  c.push_back(std::make_pair("Minor Second", minorSecond));
  c.push_back(std::make_pair("Major Second", majorSecond));
  c.push_back(std::make_pair("Minor Third", minorThird));
  c.push_back(std::make_pair("Major Third", majorThird));

  //// Triads
  c.push_back(std::make_pair("Major Triad", majorTriad));
  c.push_back(std::make_pair("Minor Triad", minorTriad));
  c.push_back(std::make_pair("Diminished Triad", diminishedTriad));
  c.push_back(std::make_pair("Augmented Triad", augmentedTriad));

  //// Suspensions
  c.push_back(std::make_pair("Sus 2", sus2));
  c.push_back(std::make_pair("Sus 4", sus4));

  //// Sevenths
  c.push_back(std::make_pair("Minor Seventh", minorSeventh));
  c.push_back(std::make_pair("Minor Major Seventh", minorMajorSeventh));
  c.push_back(std::make_pair("Dominant Seventh", dominantSeventh));
  c.push_back(std::make_pair("Major Seventh", majorSeventh));
  c.push_back(std::make_pair("Full Diminished Seventh", fullDiminishedSeventh));
  c.push_back(std::make_pair("Half Diminished Seventh", halfDiminishedSeventh));
  c.push_back(std::make_pair("Augmented Seventh", augmentedSeventh));
  c.push_back(std::make_pair("Augmented Major Seventh", augmentedMajorSeventh));

  //// Open Sevenths
  c.push_back(std::make_pair("Dominant Seventh -3", thirdlessDominantSeventh));
  c.push_back(std::make_pair("Major Seventh -3", thirdlessMajorSeventh));
  c.push_back(std::make_pair("Minor Seventh -5", fifthlessMinorSeventh));
  c.push_back(std::make_pair("Minor Major Seventh -5", fifthlessMinorMajorSeventh));
  c.push_back(std::make_pair("Dominant Seventh -5", fifthlessDominantSeventh));
  c.push_back(std::make_pair("Major Seventh -5", fifthlessMajorSeventh));

  return c;
}();
