#include <iostream>

#include <vector>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/filesystem.hpp>

#include <tinyxml2.h>

#include <sndfile.hh>

#include "SampleSorter/AbletonSample.hpp"

AbletonSample::AbletonSample(const AbletonSample & other) 
  : Sample(other)
{
  name = other.name;
  docExists = false;
  wavesExist = false;
}

AbletonSample::AbletonSample(std::string file) 
  : Sample(file) {

    name = boost::filesystem::path(file).stem().native();

    try {
      process();
    } catch (int e) {
      std::cout << "The file could not be processed" << std::endl;
    }

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

long AbletonSample::getSampleRate() {
  if (wavesExist) return sampleRate;
  getWaves();
  return sampleRate;
}

std::string AbletonSample::getName() const {
  return name;
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

  // get file path
  tinyxml2::XMLElement * pathNode = audioNode -> FirstChildElement("SampleRef")
                                      -> FirstChildElement("SourceContext")
                                      -> FirstChildElement("SourceContext")
                                      -> FirstChildElement("OriginalFileRef")
                                      -> FirstChildElement("FileRef");

  tinyxml2::XMLElement * dirNode = pathNode -> FirstChildElement("SearchHint")
                                    -> FirstChildElement("PathHint")
                                    -> FirstChildElement("RelativePathElement");

  std::string path = "/";
  while (dirNode != nullptr) {
    path += dirNode -> Attribute("Dir");
    path += "/";
    dirNode = dirNode -> NextSiblingElement("RelativePathElement");
  }

  tinyxml2::XMLElement * nameNode = pathNode -> FirstChildElement("Name");
  path += nameNode -> Attribute("Value");

  // read file
  SndfileHandle audioFile(path);
  sampleRate = audioFile.samplerate();

  if (sampleRate == 0) {
    std::cout << "The file '" << path << "' could not be read" << std::endl;
    throw 20;
  }

  // set seek to beginning
  audioFile.seek(start * sampleRate, 0);
  // number of frames to read
  long size = (end - start) * sampleRate;
  std::vector<double> rawAudioData(size * audioFile.channels());
  // read raw interleaved channels
  audioFile.read(&rawAudioData[0], size * audioFile.channels());

  // Allocate output vector
  waves.resize(audioFile.channels());
  for (int channel = 0; channel < audioFile.channels(); channel ++) {
    waves[channel].resize(size);
  }

  // De-interleave channels
  for (long i = 0; i < size * audioFile.channels(); i++) {
    waves[i % audioFile.channels()][i/audioFile.channels()] = rawAudioData[i];
  }

  wavesExist = true;
  return waves;
}
