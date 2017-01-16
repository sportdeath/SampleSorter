#include <vector>
#include <string>
#include <fstream>

#include <mad.h>

#include "SampleSorter/MPEGFile.hpp"

std::vector<std::vector<double> > MPEGFile::read(
    const std::string filePath,
    double startSeconds,
    double endSeconds) {

  std::vector<std::vector<double> > output;

  // Initialize MAD library
  struct mad_stream mad_stream;
  struct mad_frame mad_frame;
  struct mad_synth mad_synth;

  mad_stream_init(&mad_stream);
  mad_frame_init(&mad_frame);
  mad_synth_init(&mad_synth);

  // Open file
  std::ifstream file;
  file.open(filePath, std::ios::in | std::ios::binary | std::ios::ate);

  if (not file.is_open()) {
    // could not read file
    return output;
  }

  // Read the data
  unsigned long size = file.tellg();
  std::vector<unsigned char> data(size);

  file.seekg (0, std::ios::beg);
  file.read ((char *)data.data(), size);

  // data in memory, close file
  file.close();

  // Input buffer stream
  mad_stream_buffer(&mad_stream, data.data(), size);

  // Read to first frame
  while(true) {
    if(mad_frame_decode(&mad_frame, &mad_stream)) {
      if (MAD_RECOVERABLE(mad_stream.error)) {
        continue;
      } else if (mad_stream.error == MAD_ERROR_BUFLEN) {
        continue;
      } else {
        break;
      }
    } else {
      break;
    }
  }

  // Read the first frame
  mad_synth_frame(&mad_synth, &mad_frame);

  unsigned long sampleRate = mad_frame.header.samplerate;
  long numChannels = mad_synth.pcm.channels;

  // allocate output
  output.resize(numChannels);

  long startSample = startSeconds * sampleRate;
  long endSample = endSeconds * sampleRate;

  for (short channel = 0; channel < numChannels; channel++) {
    output[channel].resize(endSample - startSample);
  }

  // read the output
  long sampleOffset = 0;
  while (sampleOffset < endSample) {
    mad_synth_frame(&mad_synth, &mad_frame);

    sampleOffset = readFrame(
        &mad_synth.pcm, 
        output, 
        sampleOffset, 
        startSample, 
        endSample);

    if (mad_frame_decode(&mad_frame, &mad_stream)) {
      if (MAD_RECOVERABLE(mad_stream.error)) {
        continue;
      } else if (mad_stream.error == MAD_ERROR_BUFLEN) {
        continue;
      } else {
        break;
      }
    }
  }

  mad_synth_finish(&mad_synth);
  mad_frame_finish(&mad_frame);
  mad_stream_finish(&mad_stream);

  return output;
}

long MPEGFile::readFrame(const struct mad_pcm * pcm,
               std::vector< std::vector<double> > & output,
               long sampleOffset,
               long startSample,
               long endSample) {

  // if the last sample in the frame is after our start sample
  // and the first sample in the frame is before our end sample
  if (((sampleOffset + pcm -> length) > startSample) and
      (sampleOffset <= endSample)) {

    for (short sampleNum = 0; sampleNum < pcm -> length; sampleNum ++) {
      for (short channel = 0; channel < pcm -> channels; channel++) {

        if ((sampleNum + sampleOffset >= startSample) and
            (sampleNum + sampleOffset < endSample)) {

          output[channel][sampleNum + sampleOffset - startSample] =
            pcm -> samples[channel][sampleNum];

        }
      }
    }
  }

  return sampleOffset + pcm -> length;
}
