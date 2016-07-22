#include <string>

#include <SampleSorter/AbletonSample.hpp>

#include "gtest/gtest.h"

TEST(SamplesTest, FirstTest) {
  std::string file = "../testFile.alc";
  AbletonSample s(file);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
