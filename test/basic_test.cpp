#include <string>
#include <iostream>

#include <SampleSorter/AbletonSampleFile.hpp>

#include <boost/filesystem.hpp>


void processFiles(boost::filesystem::path pathName, std::string userLibrary) {

  // Make a list of all the samples
  std::vector<AbletonSampleFile> allSamples;
  using namespace boost::filesystem;

  // Iterate over all paths in 
  directory_iterator end_itr;
  for (directory_iterator itr(pathName);
      itr != end_itr;
      itr++) {

    if ( is_directory(itr -> status())) {

      // If it is a directory, recurse
      processFiles(itr -> path(), userLibrary);

    } else if (itr -> path().extension() == ".alc") {

      // We found a sample file!
      std::string fileName = itr -> path().native();
      std::cout << "File: " << fileName << std::endl;

      // Check the other samples to see if
      // this one is similar
      AbletonSampleFile sample(fileName, userLibrary, true);
      std::cout << sample.getReferenceFilePath() << std::endl;
      sample.process();

      int tuning = sample.getAudioSample() -> getTuningCents();
      double tempo = 60.*(sample.getAudioSample() -> getBeatWithTuning());

      std::cout << tuning << std::endl;
      std::cout << tempo << std::endl;

      //for (long i = 0; i < allSamples.size(); i++) {
      //  sample.isCompatible(allSamples[i]);
      //}

      // Add the sample to the list
      sample.writeToFile();
      allSamples.push_back(sample);
    }
  }
}

int main(int argc, char ** argv) {
  // The library file
  std::string library = "/Users/tfh/Documents/Samples/Binaural";
  // The user library file
  std::string userLibrary = "/Users/tfh/Documents/UserLibrary/";

  // Convert it to a boost path
  boost::filesystem::path libraryPath(library);

  // Process!
  processFiles(libraryPath, userLibrary);
}
