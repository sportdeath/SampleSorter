#include <vector>
#include <cmath>

class TempoFunction {
  private:
    double tempo;
    std::vector<double> aCoefficients;
    std::vector<double> bCoefficients;
    double totalSeconds;
    double theOne;
  public:
    size_t getDegree() const {
      return coefficients.size();
    }

    double getTempo(double time) const {
      double tempo = magnitudeCoefficients[0];
      for (long i = 1; i < getDegree(); i++) {
        tempo += magnitudeCoefficients[i] *
                 sin(phaseCoefficients[i] + (2*M_PI*i*time)/totalSeconds);
      }
      return tempo;
    }

    // integral from 0 to end
    double getNumBeats(double end) const {
      double numBeats = magnitudeCoefficients[0] * end;
      for (size_t i = 1; i < getDegree(); i++) {
        numBeats -= (magnitudeCoefficients[i] * totalSeconds)
                    /(2*M_PI*i) 
                    * cos(phaseCoefficients[i] + (2*M_PI*end)/totalSeconds);
      }
      return numBeats;
    }

    // integral from start to end
    double getNumBeats(double start, double end) const {
      return getNumBeats(end) - getNumBeats(start);
    }

    // the time at which the kth beat occured
    // using newtons method
    double getKthBeat(double k, double accuracy = 0.001) const {
      double guess = magnitudeCoefficients[0] * k;
      double theOneBeats = getNumBeats(theOne);

      double error = getNumBeats(guess) - theOneBeats - k;
      while (std::abs(error) > accuracy) {
        guess -= error/getTempo(guess);
        error = getNumBeats(guess) - theOneBeats - k;
      }

      return guess;
    }

    double distanceFromBeat(double time, double beatLeft, double beatRight) const {
      // this function is linear
      // f(beatLeft) = -1
      // f(beatRight) = 1
      double dist = (2./(beatLeft - beatRight)) * (time - beatRight) + 1;
      return 1-std::abs(dist);
    }

    double getValue(const std::vector<double> & onsets) const {
      double previousBeat = getKthBeat(-1); 
      double beat = theOne;
      long k = 0;

      double value = 0;
      while (previousBeat < totalSeconds) {
        // for all bins in between the two beats
        // convert the bin to seconds. 
        for (size_t bin = std::min(Units::timeToBin(previousBeat), 0);
            bin < std::max(Units::timeToBin(beat), onsets.size());
            bin ++) {
          double binSeconds = Units::binToTime(bin);
          value += onsets[bin] * distanceFromBeat(time, previousBeat, beatRight);

          k += 1;
          previousBeat = beat;
          beat = getKthBeat(k);
        }
      }
      
      return value;
    }

    
    Tempo(
        size_t degrees, 
        double guessTempo, 
        const std::vector<double> & onsets,
        double percentageError,
        double percentageStep) {

      magnitudeCoefficients.resize(degrees);
      phaseCoefficients.resize(degrees);
      magnitudeCoefficients[0] = guessTempo;

      double minValue = onsets.size();
      double minTempo = guessTempo;

      for (double trialTempo = guessTempo * (1. - percentageError);
           trialTempo < guessTempo * (1. + percentageError);
           trialTempo += percentageError * guessTempo) {

        magnitudeCoefficients[0] = trialTempo;
        double theFirstBeat = getKthBeat(1);

        for (double theOne = 0;
            theOne < Units::timeToBins(theFirstBeat);
            theOne ++) {

          double newValue = getValue(onsets);
          if (newValue < minValue) {
            minValue = newValue;
            minTempo = trialTempo;
          }
        }
      }

      for (size_t degree = 1; degree < getDegree(); degree++) {
        double minAdjustment = 0;

        for (double magnitudeGuess = - 
          coefficients[degree] = adjustment;
          if (degree == 0) {
            coefficients[degree] += guessTempo;
          }
          double theFirstBeat = getKthBeat(1);
          for (theOne = 0; 
               theOne <= Units::timeToBin(theFirstBeat);
               theOne ++) {
            double newValue = getValue(onsets);
            if (newValue < minValue) {
              minValue = newValue;
              minAdjustment = adjustment;
            }
        }
        coefficients[degree] = minAdjustment;
    }
};
