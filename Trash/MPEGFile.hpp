#ifndef MPEG_FILE_H
#define MPEG_FILE_H

#include <vector>

class MPEGFile {
  public:
    /**
     * Returns a vector whose size
     * is the number of channels. Each channel
     * contains the uncompressed audio data
     * in the range 0...1.
     */
    static std::vector<std::vector<double> > read(
        const std::string filePath,
        double startSeconds,
        double endSeoncds);

    static long readFrame(const struct mad_pcm * pcm,
                          std::vector< std::vector<double> > & ouput,
                          long sampleOffset,
                          long startSample,
                          long endSample);
};

#endif
