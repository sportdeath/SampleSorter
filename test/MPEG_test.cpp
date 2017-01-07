#include <string>
#include <vector>

#include <SampleSorter/AudioFile.hpp>
#include <Plotting/Plotting.hpp>

#include "gtest/gtest.h"

TEST(AudioFileTest, BasicRead) {
  std::string filePath = "../TestFiles/MPEGFiles/Aintnocrime.mp3";
  //std::string filePath = "../TestFiles/MPEGFiles/ComeTogether.m4a";
  //std::string filePath = "../TestFiles/MPEGFiles/WouldntItBeNice.mp3";
  //std::string filePath = "../TestFiles/MPEGFiles/ArcadeFire.mp3";
  //std::string filePath = "../TestFiles/MPEGFiles/Bear.mp3";
  //
  long sampleRate;
  
  std::vector<std::vector<double> > output = AudioFile::read(filePath, 0, 10, &sampleRate);

  Plotting::plotVector(output[0]);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
