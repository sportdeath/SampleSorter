#include <vector>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>

#include <tinyxml2.h>

#include <sndfile.hh>

#include <SampleSorter/AbletonSample.hpp>

AbletonSample::AbletonSample(std::string file) 
  : Sample(file) {

    // get rid of waves object now that
    // pre-processing is over.
    // It won't be needed again. 
    if (wavesExist) {
      waves = std::vector< std::vector<double> >();
    }
}

tinyxml2::XMLDocument * AbletonSample::getDoc() {
  if (docExists) return &doc;

  // Read file
  std::stringstream unzipped;

  std::ifstream zipped(getFile(), 
                       std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(zipped);
  boost::iostreams::copy(in, unzipped);

  // ungzip file
  doc.Parse(&unzipped.str()[0]);
  docExists = true;

  return &doc;
}


std::vector< std::vector<double> > AbletonSample::getWaves() {
  if (wavesExist) return waves;

  // get file path
  std::string path;
  // get start
  double start;
  // get end
  double end;

  // read file
  SndfileHandle audioFile(path);
  // set seek to beginning
  audioFile.seek(start * audioFile.samplerate(), 0);
  // number of frames to read
  long size = (start - end) * audioFile.samplerate();
  std::vector<double> rawAudioData(size * audioFile.channels());
  // read raw interleaved channels
  audioFile.read(&rawAudioData[0], size);

  // Allocate output vector
  std::vector< std::vector<double> > channels(audioFile.channels());
  for (int channel = 0; channel < audioFile.channels(); channel ++) {
    channels[channel].resize(size);
  }

  // De-interleave channels
  for (long i = 0; i < size * audioFile.channels(); i++) {
    channels[i % audioFile.channels()][i/audioFile.channels()] = rawAudioData[i];
  }

  wavesExist = true;
  return channels;
}
