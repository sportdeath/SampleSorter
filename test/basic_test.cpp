#include <string>
#include <iostream>

#include <SampleSorter/AbletonSample.hpp>

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"


void processFiles(boost::filesystem::path pathName) {
using namespace boost::filesystem;

  directory_iterator end_itr;
  for (directory_iterator itr(pathName);
      itr != end_itr;
      itr++) {
    if ( is_directory(itr -> status())) {
      processFiles(itr -> path());
    } else if (itr -> path().extension() == ".alc") {
      std::string fileName = itr -> path().native();
      std::cout << "Processing " << itr -> path().stem() << std::endl;
      AbletonSample sample(fileName);
    }
  }
}


TEST(SamplesTest, FirstTest) {
  std::string library = "/Users/tfh/Dropbox (MIT)/UserLibrary/SampleLibrary/Chops";
  boost::filesystem::path libraryPath(library);
  processFiles(libraryPath);
  
  //AbletonSample("../testFile.alc");
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
