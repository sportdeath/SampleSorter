#include <vector>
#include <string>

#include <mp4v2/mp4v2.h>

#include "SampleSorter/MP4File.hpp"

std::vector<std::vector<double> > MP4File::read(
    std::string filePath,
    double startSeconds,
    double endSeconds) {
  
  std::vector<std::vector<double> > output;
  MP4FileHandle file = MP4Read(filePath.c_str());

  file.

  return output;
}
