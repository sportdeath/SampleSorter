#include <iostream>

#include <vector>
#include <fstream>
#include <sstream>

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

double AbletonSampleFile::getSampleLength() const {
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

tinyxml2::XMLElement * AbletonSampleFile::getAudioNode() {
  return doc.RootElement()
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
}

tinyxml2::XMLElement * AbletonSampleFile::getLoopNode() {
  return getAudioNode() -> FirstChildElement("Loop");
}


bool AbletonSampleFile::readDoc() {

  // Start and end times in seconds
  startSeconds = getLoopNode() -> FirstChildElement("LoopStart") 
                                -> DoubleAttribute("Value");
  endSeconds = getLoopNode() -> FirstChildElement("LoopEnd") 
                                -> DoubleAttribute("Value");

  tinyxml2::XMLElement * pathNode;

  pathNode = getAudioNode() -> FirstChildElement("SampleRef")
                        -> FirstChildElement("SourceContext")
                        -> FirstChildElement("SourceContext");

  if (pathNode != nullptr) {
    pathNode = pathNode -> FirstChildElement("BrowserContentPath");

    referenceFilePath = pathNode -> Attribute("Value");
    // remove "userlibrary:"
    referenceFilePath.erase(0,11);

    // replace "%20" with spaces
    size_t replaceIndex = 0;
    while (true) {
      replaceIndex = referenceFilePath.find("%20", replaceIndex);
      if (replaceIndex == std::string::npos) break;

      referenceFilePath.replace(replaceIndex, 3, " ");
      replaceIndex += 1;
    }
    replaceIndex = 0;

    // remove the first #
    replaceIndex = referenceFilePath.find("#", replaceIndex);
    referenceFilePath.erase(replaceIndex, 1);

    // replace the :s with /
    replaceIndex = 0;
    while (true) {
      replaceIndex = referenceFilePath.find(":", replaceIndex);
      if (replaceIndex == std::string::npos) break;

      referenceFilePath.replace(replaceIndex, 1, "/");
      replaceIndex += 1;
    }
  } else { 
    pathNode = getAudioNode() -> FirstChildElement("SampleRef")
                              -> FirstChildElement("FileRef");

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
  }

  // Check for preprocessing
  tinyxml2::XMLElement * sortData = doc.RootElement() 
                        -> FirstChildElement("LiveSet")
                        -> FirstChildElement("SampleSorter");

  if (sortData != nullptr) {
    long tuningCents = getAudioNode() -> FirstChildElement("PitchFine") 
                                      -> IntAttribute("Value");
    double rawBeat = sortData -> DoubleAttribute("RawBeat");
    double theOne = sortData -> DoubleAttribute("TheOne");

    std::vector<Octave> chords;

    tinyxml2::XMLElement * chordElement = sortData -> FirstChildElement("Chord");
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

    sample = AudioSample(tuningCents, rawBeat, theOne, endSeconds - startSeconds, chords);

    return true;
  }

  return false;
}

void AbletonSampleFile::writeToFile() {

  tinyxml2::XMLElement * sortData;

  sortData = doc.RootElement() -> FirstChildElement("LiveSet")
                               -> FirstChildElement("SampleSorter");

  doc.RootElement() -> FirstChildElement("LiveSet")
                    -> DeleteChild(sortData);

  sortData = doc.NewElement("SampleSorter");

  sortData -> SetAttribute("ReferenceFilePath", referenceFilePath.c_str());

  sortData -> SetAttribute("RawBeat", getAudioSample() -> getBeatRaw());

  getLoopNode() -> FirstChildElement("HiddenLoopStart") 
    -> SetAttribute("Value", startSeconds);

  getLoopNode() -> FirstChildElement("HiddenLoopEnd") 
    -> SetAttribute("Value", endSeconds);

  getAudioNode() -> FirstChildElement("Name") -> DeleteAttribute("Name");
  getAudioNode() -> FirstChildElement("Name")
    -> SetAttribute("Value", (getFileName() + " @" 
        + std::to_string(60. * (getAudioSample() -> getBeatWithTuning()))).c_str());

  sortData -> SetAttribute("TheOne", getAudioSample() -> getTheOneRaw());

  int tuningCents = getAudioSample() -> getTuningCents();
  getAudioNode() -> FirstChildElement("PitchFine")
                 -> SetAttribute("Value", tuningCents);

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

  doc.RootElement() -> FirstChildElement("LiveSet") -> InsertEndChild(sortData); 

  // make a string
  tinyxml2::XMLPrinter printer;
  doc.Accept( &printer);

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
