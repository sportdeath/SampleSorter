#include <string>
#include <iostream>

#include <SampleSorter/AbletonSampleFile.hpp>

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"


void processFiles(boost::filesystem::path pathName) {

  std::vector<AbletonSampleFile> allSamples;
using namespace boost::filesystem;

  directory_iterator end_itr;
  for (directory_iterator itr(pathName);
      itr != end_itr;
      itr++) {
    if ( is_directory(itr -> status())) {
      processFiles(itr -> path());
    } else if (itr -> path().extension() == ".alc") {
      std::string fileName = itr -> path().native();
      std::cout << allSamples.size() << std::endl;
      AbletonSampleFile sample(fileName);
      for (long i = 0; i < allSamples.size(); i++) {
        //sample.isCompatible(allSamples[i]);
      }
      allSamples.push_back(sample);
      std::cout << "new size: " << allSamples.size() << std::endl;
    }
  }
}


TEST(SamplesTest, FirstTest) {
  //std::string library = "/Users/tfh/Dropbox (MIT)/UserLibrary/SampleLibrary/Chops";
  //boost::filesystem::path libraryPath(library);
  //processFiles(libraryPath);
  
  AbletonSampleFile a("../testFile.alc");
  a.process();
  std::cout << "tuning cents: " << a.getAudioSample() -> getTuningCents() << std::endl;
  std::cout << "tempo: " << 60.*(a.getAudioSample() -> getBeatWithTuning()) << std::endl;
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
