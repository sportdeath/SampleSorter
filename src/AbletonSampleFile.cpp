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

#include "SampleSorter/SampleFile.hpp"
#include "SampleSorter/AbletonSampleFile.hpp"

AbletonSampleFile::AbletonSampleFile(const AbletonSampleFile & other) 
  : SampleFile(other)
{
  AbletonSampleFile(other.filePath);
}

AbletonSampleFile::AbletonSampleFile(std::string filePath) 
  : SampleFile(filePath) {
    getDoc();
    readDoc();
}

std::string AbletonSampleFile::getReferenceFilePath() const {
  return referenceFilePath;
}

std::string AbletonSampleFile::getReferenceFileName() const {
  return boost::filesystem::path(referenceFilePath).stem().native();
}

double AbletonSampleFile::getSampleLength() const {
  return endSeconds - startSeconds;
}

long AbletonSampleFile::getSampleRate() const {
  return audioFile.samplerate();
}

std::vector< std::vector<double> > AbletonSampleFile::getWaves() {
  // set seek to beginning
  audioFile.seek(startSeconds * getSampleRate(), 0);
  // number of frames to read
  long size = (endSeconds - startSeconds) * getSampleRate();
  std::vector<double> rawAudioData(size * audioFile.channels());
  // read raw interleaved channels
  audioFile.read(&rawAudioData[0], size * audioFile.channels());

  // Allocate output vector
  std::vector<std::vector<double> > waves(audioFile.channels());
  for (int channel = 0; channel < audioFile.channels(); channel ++) {
    waves[channel].resize(size);
  }

  // De-interleave channels
  for (long i = 0; i < size * audioFile.channels(); i++) {
    waves[i % audioFile.channels()][i/audioFile.channels()] = rawAudioData[i];
  }

  return waves;
}

void AbletonSampleFile::getDoc() {
  // Ungzip
  std::stringstream unzipped;

  std::ifstream zipped(filePath, 
                       std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(zipped);
  boost::iostreams::copy(in, unzipped);

  // convert from XML
  doc.Parse(&unzipped.str()[0]);
}

void AbletonSampleFile::readDoc() {
  tinyxml2::XMLElement * audioNode = doc.RootElement()
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
  startSeconds = loopNode -> FirstChildElement("LoopStart") 
                                -> DoubleAttribute("Value");
  endSeconds = loopNode -> FirstChildElement("LoopEnd") 
                                -> DoubleAttribute("Value");

  tinyxml2::XMLElement * pathNode = audioNode -> FirstChildElement("SampleRef")
                                              -> FirstChildElement("FileRef");
  if (pathNode == nullptr) {
    pathNode = audioNode -> FirstChildElement("SampleRef")
                         -> FirstChildElement("SourceContext")
                         -> FirstChildElement("SourceContext")
                         -> FirstChildElement("OriginalFileRef")
                         -> FirstChildElement("FileRef");
  }

  tinyxml2::XMLElement * dirNode = pathNode -> FirstChildElement("SearchHint")
                                    -> FirstChildElement("PathHint")
                                    -> FirstChildElement("RelativePathElement");

  // get file path
  referenceFilePath = "/";
  while (dirNode != nullptr) {
    referenceFilePath += dirNode -> Attribute("Dir");
    referenceFilePath += "/";
    dirNode = dirNode -> NextSiblingElement("RelativePathElement");
  }

  // get file name
  tinyxml2::XMLElement * nameNode = pathNode -> FirstChildElement("Name");
  referenceFilePath += nameNode -> Attribute("Value");

  audioFile = SndfileHandle(referenceFilePath);
}
