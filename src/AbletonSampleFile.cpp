#include <iostream>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <initializer_list>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/filesystem.hpp>

#include <tinyxml2.h>

#include "SampleSorter/SampleFile.hpp"
#include "SampleSorter/AbletonSampleFile.hpp"
#include "SampleSorter/ProcessingException.hpp"
#include "SampleSorter/AudioFile.hpp"

AbletonSampleFile::AbletonSampleFile(const AbletonSampleFile & other) 
  : SampleFile(other)
{
  AbletonSampleFile(other.filePath);
}

AbletonSampleFile::AbletonSampleFile(std::string filePath) 
  : SampleFile(filePath) {
}

bool AbletonSampleFile::readMetaData() {
  getDoc();
  return readDoc();
}

std::string AbletonSampleFile::getReferenceFilePath() const {
  return referenceFilePath;
}

std::string AbletonSampleFile::getReferenceFileName() const {
  return boost::filesystem::path(referenceFilePath).stem().native();
}

double AbletonSampleFile::getSampleSeconds() const {
  return endSeconds - startSeconds;
}

std::vector< std::vector<double> > AbletonSampleFile::extractAudio(long * sampleRate) {
  return AudioFile::read(referenceFilePath, startSeconds, endSeconds, sampleRate);
}

void AbletonSampleFile::getDoc() {
  std::stringstream unzipped;
  try {
    // Ungzip
    std::ifstream zipped(filePath, 
                         std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(zipped);
    boost::iostreams::copy(in, unzipped);
  } catch (boost::iostreams::gzip_error & e) {
    throw ProcessingException("Could not ungzip Ableton file!");
  }

  // convert from XML
  doc.Parse(&unzipped.str()[0]);
}

void AbletonSampleFile::setDoc() {
  // make a string
  tinyxml2::XMLPrinter printer;
  doc.Accept(&printer);

  // gzip
  try {
    std::stringstream unzipped;
    unzipped << printer.CStr();

    std::ofstream zipped(filePath, 
                         std::ios_base::out | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> filter;
    filter.push(boost::iostreams::gzip_compressor());
    filter.push(unzipped);
    boost::iostreams::copy(filter, zipped);
  } catch (boost::iostreams::gzip_error & e) {
    throw ProcessingException("Could not gzip Ableton file!");
  }
}

tinyxml2::XMLElement * AbletonSampleFile::getNode(
    std::initializer_list<std::string> nodes
    ) {
  tinyxml2::XMLElement * node = doc.RootElement();
  for (auto it = nodes.begin(); it < nodes.end(); it++) {
    node = node -> FirstChildElement((*it).c_str());
  }
  return node;
}

tinyxml2::XMLElement * AbletonSampleFile::getAudioNode() {
  return getNode({
      "LiveSet",
      "Tracks",
      "AudioTrack",
      "DeviceChain",
      "MainSequencer",
      "ClipSlotList",
      "ClipSlot",
      "ClipSlot",
      "Value",
      "AudioClip"
      });
}

tinyxml2::XMLElement * AbletonSampleFile::getLoopNode() {
  return getAudioNode() -> FirstChildElement("Loop");
}

tinyxml2::XMLElement * AbletonSampleFile::getSortDataNode() {
  return getNode({"LiveSet","SampleSorter"});
}

tinyxml2::XMLElement * AbletonSampleFile::getLoopStartNode() {
  return getLoopNode() -> FirstChildElement("LoopStart");
}

tinyxml2::XMLElement * AbletonSampleFile::getLoopEndNode() {
  return getLoopNode() -> FirstChildElement("LoopEnd");
}

tinyxml2::XMLElement * AbletonSampleFile::getHiddenLoopStartNode() {
  return getLoopNode() -> FirstChildElement("HiddenLoopStart");
}

tinyxml2::XMLElement * AbletonSampleFile::getHiddenLoopEndNode() {
  return getLoopNode() -> FirstChildElement("HiddenLoopEnd");
}

tinyxml2::XMLElement * AbletonSampleFile::getNameNode() {
  return getAudioNode() -> FirstChildElement("Name");
}

tinyxml2::XMLElement * AbletonSampleFile::getPitchFineNode() {
  return getAudioNode() -> FirstChildElement("PitchFine");
}

std::string AbletonSampleFile::fetchReferenceFilePath() {
  // Some awful reverse engineering
  std::string path;

  tinyxml2::XMLElement * pathNode = getAudioNode() 
    -> FirstChildElement("SampleRef")
    -> FirstChildElement("SourceContext")
    -> FirstChildElement("SourceContext");

  if (pathNode != nullptr) {
    pathNode = pathNode -> FirstChildElement("BrowserContentPath");

    path = pathNode -> Attribute("Value");
    // remove "userlibrary:"
    path.erase(0,11);

    // replace "%20" with spaces
    size_t replaceIndex = 0;
    while (true) {
      replaceIndex = path.find("%20", replaceIndex);
      if (replaceIndex == std::string::npos) break;

      path.replace(replaceIndex, 3, " ");
      replaceIndex += 1;
    }
    replaceIndex = 0;

    // remove the first #
    replaceIndex = path.find("#", replaceIndex);
    path.erase(replaceIndex, 1);

    // replace the :s with /
    replaceIndex = 0;
    while (true) {
      replaceIndex = path.find(":", replaceIndex);
      if (replaceIndex == std::string::npos) break;

      path.replace(replaceIndex, 1, "/");
      replaceIndex += 1;
    }
  } else { 
    pathNode = getAudioNode() -> FirstChildElement("SampleRef")
                              -> FirstChildElement("FileRef");

    tinyxml2::XMLElement * dirNode = pathNode -> FirstChildElement("SearchHint")
                                      -> FirstChildElement("PathHint")
                                      -> FirstChildElement("RelativePathElement");

    // get file path
    path = "/";
    while (dirNode != nullptr) {
      path += dirNode -> Attribute("Dir");
      path += "/";
      dirNode = dirNode -> NextSiblingElement("RelativePathElement");
    }

    // get file name
    tinyxml2::XMLElement * nameNode = pathNode -> FirstChildElement("Name");
    path += nameNode -> Attribute("Value");
  }

  return path;
}

bool AbletonSampleFile::readDoc() {

  // Start and end times in seconds
  startSeconds = getLoopStartNode() -> DoubleAttribute("Value");
  endSeconds = getLoopEndNode() -> DoubleAttribute("Value");

  // if this has already been processed
  if (getSortDataNode() != nullptr) {
    // Path
    referenceFilePath = getSortDataNode() -> Attribute("ReferenceFilePath");
    // Tuning
    long tuningCents = getPitchFineNode() -> IntAttribute("Value");
    // Tempo
    double rawBeat = getSortDataNode() -> DoubleAttribute("RawBeat");
    // One
    double theOne = getSortDataNode() -> DoubleAttribute("TheOne");
    // SampleRate
    double sampleRate = getSortDataNode() -> IntAttribute("SampleRate");

    // Chords
    std::vector<Octave> chords;
    tinyxml2::XMLElement * chordElement = getSortDataNode() 
      -> FirstChildElement("Chord");
    while (chordElement != nullptr) {
      std::vector<double> spectrogram(12);
      tinyxml2::XMLElement * noteElement = 
        chordElement -> FirstChildElement("Note");
      for (long i = 0; i < 12; i++) {
        spectrogram[i] = noteElement -> DoubleAttribute("Value");
        noteElement = noteElement -> NextSiblingElement("Note");
      }
      chordElement = chordElement -> NextSiblingElement("Chord");
      Octave newChord(spectrogram);
      chords.push_back(newChord);
    }

    // Make the sample!
    sample = AudioSample(
        tuningCents, 
        rawBeat, 
        theOne, 
        endSeconds - startSeconds, 
        sampleRate,
        chords
        );

    return true;
  } else {
    // it was not already processed, return false
    referenceFilePath = fetchReferenceFilePath();
    return false;
  }
}

void AbletonSampleFile::writeToFile() {
  // Delete old data to write new
  tinyxml2::XMLElement * sortData = getSortDataNode();
  doc.RootElement() -> FirstChildElement("LiveSet")
                    -> DeleteChild(sortData);
  sortData = doc.NewElement("SampleSorter");

  // Path
  sortData -> SetAttribute("ReferenceFilePath", referenceFilePath.c_str());
  // Tempo
  sortData -> SetAttribute("RawBeat", getAudioSample() -> getBeatRaw());
  // The one
  sortData -> SetAttribute("TheOne", getAudioSample() -> getTheOneRaw());
  // Sample rate
  sortData -> SetAttribute("SampleRate", int(getAudioSample() -> getSampleRate()));

  // Loop ends
  getHiddenLoopStartNode() -> SetAttribute("Value", startSeconds);
  getHiddenLoopEndNode() -> SetAttribute("Value", endSeconds);

  // Tuning
  getPitchFineNode() 
    -> SetAttribute("Value", int(getAudioSample() -> getTuningCents()));

  // Name
  getNameNode() -> SetAttribute("Value",
      (
        getFileName() + " @" +
        std::to_string(60. * (getAudioSample() -> getBeatWithTuning()))
      ).c_str()
    );

  // Chords
  std::vector<Octave> chords = getAudioSample() -> getChords();
  for (long i = 0; i < chords.size(); i++) {
    tinyxml2::XMLElement * chordElement = doc.NewElement("Chord");
    for (short j = 0; j < 12; j++) {
      tinyxml2::XMLElement * noteElement = doc.NewElement("Note");
      noteElement -> SetAttribute("Value", chords[i].getSpectrogram()[j]);
      chordElement -> InsertEndChild(noteElement);
    }
    sortData -> InsertEndChild(chordElement);
  }

  // Put in doc
  doc.RootElement() -> FirstChildElement("LiveSet") -> InsertEndChild(sortData); 
  setDoc();
}
