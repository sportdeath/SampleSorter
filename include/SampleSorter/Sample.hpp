#include <string>

class Sample {
  private:
    std::string file;

    bool isHarmonic_ = false;
    //std::vector<Chord> chords;

    bool isTuned = false;

    bool hasBeat_ = false;
    double beat;

    /**
     * Saves any processing done back to
     * the file.
     */
    //virtual void saveToFile() = 0;

    /**
     * Reads the samples values from a file
     * The result is a vector of channels,
     * where each channel is a vector or
     * sample values.
     */
    virtual std::vector< std::vector<double> > getWaves() = 0;

    virtual long getSampleRate() = 0;

    /**
     * Tunes the sample to A = 440
     * and updates the file structure
     */
    void tune();
    void findBeat();
    //void findChords();
  protected:
    void process();
  public:
    /**
     * Reads the sample information
     * from a file. If necessary, some processing is
     * done:
     *  The sample is tuned
     *  The beat is found
     *  The chords are generated
     * This result is saved for future use.
     */
    Sample(std::string file);

    std::string getFile();

    /**
     * Returns whether or not the sample
     * is harmonic or inharmonic.
     */
    bool isHarmonic();

    /**
     * Returns whether or not the sample
     * has a beat
     */
    bool hasBeat();

    /**
     * If the sample is inharmonic, this is the number
     * of rhythmically similar structures per second as
     * defined by beat onsets (typically measures).
     * Otherwise this is the minimum of the that value
     * and the maximum cyclic beat that divides the sample
     * into harmonically similar chunks ("Chords") - an
     * integer fraction of the other number.
     *
     * The beat is in units of beats per second (BPS)
     *
     * An error is returned if the sample has no beat
     */
    double getBeat();

    /**
     * Divides the sample into harmonically similar
     * chunks of equal size according the BPS and
     * returns them as Chords which contain pitch 
     * information.
     *
     * If the sample is inharmonic an error is returned
     */
    //std::vector<Chord> getChords();
};
