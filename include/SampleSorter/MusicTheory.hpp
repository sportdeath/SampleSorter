#ifndef MUSIC_THEORY_H
#define MUSIC_THEORY_H

#include "SampleSorter/Octave.hpp"

class MusicTheory {
  public:
    static const std::vector<std::pair<std::string, Octave> > chords;

    static const std::vector<std::string> noteNames;

    // Intervals
    static const Octave unison;
    static const Octave minorSecond;
    static const Octave majorSecond;
    static const Octave minorThird;
    static const Octave majorThird;
    static const Octave perfectFifth;
    
    // Triads
    static const Octave majorTriad;
    static const Octave augmentedTriad;
    static const Octave minorTriad;
    static const Octave diminishedTriad;

    // Suspensions
    static const Octave sus2;
    static const Octave sus4;

    // Sevenths
    static const Octave fullDiminishedSeventh;
    static const Octave halfDiminishedSeventh;
    static const Octave minorSeventh;
    static const Octave minorMajorSeventh;
    static const Octave dominantSeventh;
    static const Octave majorSeventh;
    static const Octave augmentedSeventh;
    static const Octave augmentedMajorSeventh;

    // Open Sevenths
    static const Octave thirdlessDominantSeventh;
    static const Octave thirdlessMajorSeventh;
    static const Octave fifthlessMinorSeventh;
    static const Octave fifthlessMinorMajorSeventh;
    static const Octave fifthlessDominantSeventh;
    static const Octave fifthlessMajorSeventh;
};

#endif
