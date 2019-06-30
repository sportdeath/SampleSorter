#include <iostream>
#include <ciso646>

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
  AbletonSampleFile(other.filePath, other.userLibrary, other.wasOrWillBeProcessed);
}

AbletonSampleFile::AbletonSampleFile(
    std::string filePath,
    std::string userLibrary,
    bool forceReprocess_) 
  : SampleFile(filePath),
    userLibrary(userLibrary)
{
  forceReprocess = forceReprocess_;
}

bool AbletonSampleFile::readMetaData() {
  getDoc();
  wasOrWillBeProcessed = (getSortDataNode() == nullptr) or forceReprocess;
  if (not wasOrWillBeProcessed)  {
    wasOrWillBeProcessed = wasOrWillBeProcessed or 
    ((getSortDataNode() -> FirstChildElement("Octave")) == nullptr);
  }
  return readDoc();
}

std::string AbletonSampleFile::getReferenceFilePath() const {
  return referenceFilePath;
}

std::string AbletonSampleFile::getReferenceFileName() const {
  return boost::filesystem::path(referenceFilePath).stem().make_preferred().string();
	  // native();
}

double AbletonSampleFile::getSampleSeconds() const {
  return endSeconds - startSeconds;
}

std::vector< std::vector<double> > AbletonSampleFile::extractAudio(long * sampleRate) {
  return AudioFile::read(referenceFilePath, startSeconds, endSeconds, sampleRate);
}

void AbletonSampleFile::getDoc() {
  std::stringstream unzipped;
  std::string xml;
  try {
    // Ungzip
    std::ifstream zipped(filePath, 
                         std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(zipped);
    boost::iostreams::copy(in, unzipped);
    xml = unzipped.str();
  } catch (boost::iostreams::gzip_error & e) {
    // Try to open it as text
    std::ifstream text(filePath);
    if (text.good()) {
      std::string firstLine;
      std::getline(text, firstLine);
      if (firstLine == "<?xml version=\"1.0\" encoding=\"UTF-8\"?>") {
        // Go back to the beginning and read the whole file
        text.clear();
        text.seekg(0, std::ios::beg);
        xml = std::string((std::istreambuf_iterator<char>(text)),
             std::istreambuf_iterator<char>());
      }
    }
    text.close();
  }

  if (xml == "") {
    throw ProcessingException("Could not read Ableton file!");
  }

  // convert from XML
  doc.Parse(&xml[0]);
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
    std::cout << std::endl;
    std::cerr << e.what() << std::endl;
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
  
  // Build the relative path
  tinyxml2::XMLElement * pathNode = 
    getAudioNode() -> FirstChildElement("SampleRef")
    -> FirstChildElement("FileRef");

  tinyxml2::XMLElement * dirNode = 
    pathNode 
    -> FirstChildElement("RelativePath")
    -> FirstChildElement("RelativePathElement");

  // get file path
  std::string path = userLibrary;
  while (dirNode != nullptr) {
    path += dirNode -> Attribute("Dir");
    path += "/";
    dirNode = dirNode -> NextSiblingElement("RelativePathElement");
  }

  // get file name
  tinyxml2::XMLElement * nameNode = 
    pathNode -> FirstChildElement("Name");
  path += nameNode -> Attribute("Value");
  // }

  return path;
}

bool AbletonSampleFile::readDoc() {
  // Start and end times in seconds
  startSeconds = getLoopStartNode() -> DoubleAttribute("Value");
  endSeconds = getLoopEndNode() -> DoubleAttribute("Value");

  referenceFilePath = fetchReferenceFilePath();

  // if this has already been processed
  if (not wasOrWillBeProcessed) {
    // Tuning
    long tuningCents = getPitchFineNode() -> IntAttribute("Value");
    // Fundemental
    short fundemental = getSortDataNode() -> IntAttribute("Fundemental");
    // Tempo
    double rawBeat = getSortDataNode() -> DoubleAttribute("RawBeat");
    // One
    double theOne = getSortDataNode() -> DoubleAttribute("TheOne");
    // SampleRate
    double sampleRate = getSortDataNode() -> IntAttribute("SampleRate");

    // Octave
    tinyxml2::XMLElement * octaveElement = getSortDataNode() 
      -> FirstChildElement("Octave");
    std::vector<double> octaveSpectrogram(12);
    for (int i = 0; i < 12; i++) {
      octaveSpectrogram[i] = octaveElement -> DoubleAttribute("Value");
      octaveElement = octaveElement -> NextSiblingElement("Octave");
    }
    Octave octave(octaveSpectrogram);

    // Chords
    std::vector<Octave> chords;
    tinyxml2::XMLElement * chordElement = getSortDataNode() 
      -> FirstChildElement("Chord");
    while (chordElement != nullptr) {
      std::vector<double> spectrogram(12);
      tinyxml2::XMLElement * noteElement = 
        chordElement -> FirstChildElement("Note");
      for (int i = 0; i < 12; i++) {
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
        fundemental,
        rawBeat, 
        theOne, 
        endSeconds - startSeconds, 
        sampleRate,
        chords,
        octave
        );

    return true;
  } else {
      // it was not already processed, return false
    return false;
  }
}

void AbletonSampleFile::writeToFile() {
  if (not wasOrWillBeProcessed) return;

  // Delete old data if it exists
  tinyxml2::XMLElement * sortData = getSortDataNode();
  if (sortData != NULL) {
    doc.RootElement() -> FirstChildElement("LiveSet")
                      -> DeleteChild(sortData);
  }
  // Make a new node
  sortData = doc.NewElement("SampleSorter");

  // Path
  sortData -> SetAttribute("ReferenceFilePath", referenceFilePath.c_str());
  // Tempo
  sortData -> SetAttribute("RawBeat", getAudioSample() -> getBeatRaw());
  // The one
  sortData -> SetAttribute("TheOne", getAudioSample() -> getTheOneRaw());
  // Sample rate
  sortData -> SetAttribute("SampleRate", int(getAudioSample() -> getSampleRate()));
  // Fundemental
  sortData -> SetAttribute("Fundemental", getAudioSample() -> getFundemental());

  // Loop ends
  // set the loop start to be the first beat
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

  // Octave
  Octave octave = getAudioSample() -> getOctave();
  for (short i = 0; i < 12; i++) {
    tinyxml2::XMLElement * octaveElement = doc.NewElement("Octave");
    octaveElement -> SetAttribute("Value", octave.getSpectrogram()[i]);
    sortData -> InsertEndChild(octaveElement);
  }

  // Chords
  std::vector<Octave> chords = getAudioSample() -> getChords();
  for (size_t i = 0; i < chords.size(); i++) {
    tinyxml2::XMLElement * chordElement = doc.NewElement("Chord");
    for (int j = 0; j < 12; j++) {
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
