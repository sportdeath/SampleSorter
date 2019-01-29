#include <iostream>
#include <string>

#include <sample_sorter/audio.hpp>

int main() {
  std::string filename = "/home/tfh/sample_sorter/test/smokey.wav";
  double sample_rate;

  std::vector<std::vector<double>> audio =
    sample_sorter::Audio::read(filename, &sample_rate);

  std::cout << "Sample rate: " << sample_rate << std::endl;
  std::cout << "Number of channels: " << audio.size() << std::endl;
  std::cout << "Number of samples: " << audio[0].size() << std::endl;
  std::cout << "Audio length: " << audio[0].size()/sample_rate << std::endl;
}
