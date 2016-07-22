#include <iostream>

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

    process();

    // get rid of waves object now that
    // pre-processing is over.
    // It won't be needed again. 
    if (wavesExist) {
      waves = std::vector< std::vector<double> >();
    }
}

tinyxml2::XMLDocument * AbletonSample::getDoc() {
  if (docExists) return &doc;

  // Ungzip
  std::stringstream unzipped;

  std::ifstream zipped(getFile(), 
                       std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(zipped);
  boost::iostreams::copy(in, unzipped);

  // convert from xml
  doc.Parse(&unzipped.str()[0]);
  docExists = true;

  return &doc;
}


std::vector< std::vector<double> > AbletonSample::getWaves() {
  if (wavesExist) return waves;

  tinyxml2::XMLElement * audioNode = getDoc() -> RootElement()
                                          -> FirstChildElement("LiveSet")
                                          -> FirstChildElement("Tracks")
                                          -> FirstChildElement("AudioTrack")
                                          -> FirstChildElement("DeviceChain")
                                          -> FirstChildElement("MainSequencer")
                                          -> FirstChildElement("ClipSlotList")
                                          -> FirstChildElement("ClipSlot")
                                          -> FirstChildElement("ClipSlot")
                                          -> FirstChildElement("Value")
                                          -> FirstChildElement("AudioClip");

  tinyxml2::XMLElement * loopNode = audioNode -> FirstChildElement("Loop");


  // Start and end times in seconds
  double start = loopNode -> FirstChildElement("LoopStart") 
                                -> DoubleAttribute("Value");
  double end = loopNode -> FirstChildElement("LoopEnd") 
                                -> DoubleAttribute("Value");

  tinyxml2::XMLElement * pathNode = audioNode -> FirstChildElement("SampleRef")
                                      -> FirstChildElement("SourceContext")
                                      -> FirstChildElement("SourceContext")
                                      -> FirstChildElement("OriginalFileRef")
                                      -> FirstChildElement("FileRef");

  tinyxml2::XMLElement * dirNode = pathNode -> FirstChildElement("SearchHint")
                                    -> FirstChildElement("PathHint");

  tinyxml2::XMLElement * pathList = dirNode -> FirstChildElement("RelativePathElement");

  std::string path = "/";
  while (pathList != nullptr) {
    path += pathList -> Attribute("Dir");
    path += "/";
    pathList = pathList -> NextSiblingElement("RelativePathElement");
  }

  tinyxml2::XMLElement * nameNode = pathNode -> FirstChildElement("Name");
  path += nameNode -> Attribute("Value");

  // get file path

  // read file
  SndfileHandle audioFile(path);
  // set seek to beginning
  audioFile.seek(start * audioFile.samplerate(), 0);
  // number of frames to read
  long size = (end - start) * audioFile.samplerate();
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
