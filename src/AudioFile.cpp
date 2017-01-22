diff --git a/.gitignore b/.gitignore
index 197eb5c..a9d338f 100644
--- a/.gitignore
+++ b/.gitignore
@@ -1,5 +1,8 @@
 # The Cmake build directory
 build/
+TestFiles/
+TestFilesBackup.zip
+python_src/
 
 # OS generated files #
 ######################
diff --git a/CMakeLists.txt b/CMakeLists.txt
index 9b2526c..6e038ff 100755
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -19,6 +19,8 @@ set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/modules/)
 set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
 set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
 set(PROJECT_SOURCES
+  ${PROJECT_SOURCE_DIR}/AudioFile.cpp
+
   ${PROJECT_SOURCE_DIR}/AbletonSampleFile.cpp
   ${PROJECT_SOURCE_DIR}/SampleFile.cpp
 
@@ -28,7 +30,9 @@ set(PROJECT_SOURCES
   ${PROJECT_SOURCE_DIR}/Octave.cpp
   ${PROJECT_SOURCE_DIR}/Tempo.cpp
   ${PROJECT_SOURCE_DIR}/EqualLoudness.cpp
+
   ${PROJECT_SOURCE_DIR}/MusicTheory.cpp
+  ${PROJECT_SOURCE_DIR}/Units.cpp
 
   ${PROJECT_SOURCE_DIR}/PythonExternals.cpp
 
@@ -66,6 +70,27 @@ if(NOT LIBSNDFILE_FOUND)
   message(ERROR " LIBSNDFILE not found!")
 endif(NOT LIBSNDFILE_FOUND)
 
+# LibMad
+find_package(LIBMAD REQUIRED)
+include_directories(${LIBMAD_INCLUDE_DIR})
+set(LIBS ${LIBS} ${LIBMAD_LIBRARY})
+if(NOT LIBMAD_FOUND)
+  message(ERROR "LIBMAD not found!")
+endif(NOT LIBMAD_FOUND)
+
+# FFMPEG
+find_package(FFMPEG REQUIRED)
+include_directories(${FFMPEG_INCLUDE_DIR})
+set(LIBS 
+  ${LIBS} 
+  ${FFMPEG_LIBRARIES} 
+  #${FFMPEG_avformat_LIBRARY}
+  #${FFMPEG_avcodec_LIBRARY}
+)
+if(NOT FFMPEG_FOUND)
+  message(ERROR "FFMPEG not found!")
+endif(NOT FFMPEG_FOUND)
+
 # FFTW3
 find_package(FFTW REQUIRED)
 include_directories(${FFTW_INCLUDES})
@@ -74,7 +99,6 @@ if(NOT FFTW_FOUND)
   message(ERROR " FFTW not found!")
 endif(NOT FFTW_FOUND)
 
-
 ####################
 ## Library Creation
 ####################
@@ -117,3 +141,5 @@ macro(MAKE_TEST test_name)
 endmacro()
 
 MAKE_TEST(basic_test)
+MAKE_TEST(MPEG_test)
+MAKE_TEST(tensorFlow_test)
diff --git a/ComeTogether.m4a.reapeaks b/ComeTogether.m4a.reapeaks
new file mode 100644
index 0000000..6d61fd6
Binary files /dev/null and b/ComeTogether.m4a.reapeaks differ
diff --git a/Trash/MP4File.cpp b/Trash/MP4File.cpp
new file mode 100644
index 0000000..5df658d
--- /dev/null
+++ b/Trash/MP4File.cpp
@@ -0,0 +1,19 @@
+#include <vector>
+#include <string>
+
+#include <mp4v2/mp4v2.h>
+
+#include "SampleSorter/MP4File.hpp"
+
+std::vector<std::vector<double> > MP4File::read(
+    std::string filePath,
+    double startSeconds,
+    double endSeconds) {
+  
+  std::vector<std::vector<double> > output;
+  MP4FileHandle file = MP4Read(filePath.c_str());
+
+  file.
+
+  return output;
+}
diff --git a/Trash/MP4File.hpp b/Trash/MP4File.hpp
new file mode 100644
index 0000000..f124804
--- /dev/null
+++ b/Trash/MP4File.hpp
@@ -0,0 +1,12 @@
+#ifndef MP4_FILE_H
+#define MP4_FILE_H
+
+class MP4File {
+  public:
+    static std::vector<std::vector<double> >
+      read(std::string filePath,
+           double startSeconds,
+           double endSeconds);
+};
+
+#endif
diff --git a/Trash/MPEGFile.cpp b/Trash/MPEGFile.cpp
new file mode 100644
index 0000000..a8651e0
--- /dev/null
+++ b/Trash/MPEGFile.cpp
@@ -0,0 +1,134 @@
+#include <vector>
+#include <string>
+#include <fstream>
+
+#include <mad.h>
+
+#include "SampleSorter/MPEGFile.hpp"
+
+std::vector<std::vector<double> > MPEGFile::read(
+    const std::string filePath,
+    double startSeconds,
+    double endSeconds) {
+
+  std::vector<std::vector<double> > output;
+
+  // Initialize MAD library
+  struct mad_stream mad_stream;
+  struct mad_frame mad_frame;
+  struct mad_synth mad_synth;
+
+  mad_stream_init(&mad_stream);
+  mad_frame_init(&mad_frame);
+  mad_synth_init(&mad_synth);
+
+  // Open file
+  std::ifstream file;
+  file.open(filePath, std::ios::in | std::ios::binary | std::ios::ate);
+
+  if (not file.is_open()) {
+    // could not read file
+    return output;
+  }
+
+  // Read the data
+  unsigned long size = file.tellg();
+  std::vector<unsigned char> data(size);
+
+  file.seekg (0, std::ios::beg);
+  file.read ((char *)data.data(), size);
+
+  // data in memory, close file
+  file.close();
+
+  // Input buffer stream
+  mad_stream_buffer(&mad_stream, data.data(), size);
+
+  // Read to first frame
+  while(true) {
+    if(mad_frame_decode(&mad_frame, &mad_stream)) {
+      if (MAD_RECOVERABLE(mad_stream.error)) {
+        continue;
+      } else if (mad_stream.error == MAD_ERROR_BUFLEN) {
+        continue;
+      } else {
+        break;
+      }
+    } else {
+      break;
+    }
+  }
+
+  // Read the first frame
+  mad_synth_frame(&mad_synth, &mad_frame);
+
+  unsigned long sampleRate = mad_frame.header.samplerate;
+  long numChannels = mad_synth.pcm.channels;
+
+  // allocate output
+  output.resize(numChannels);
+
+  long startSample = startSeconds * sampleRate;
+  long endSample = endSeconds * sampleRate;
+
+  for (short channel = 0; channel < numChannels; channel++) {
+    output[channel].resize(endSample - startSample);
+  }
+
+  // read the output
+  long sampleOffset = 0;
+  while (sampleOffset < endSample) {
+    mad_synth_frame(&mad_synth, &mad_frame);
+
+    sampleOffset = readFrame(
+        &mad_synth.pcm, 
+        output, 
+        sampleOffset, 
+        startSample, 
+        endSample);
+
+    if (mad_frame_decode(&mad_frame, &mad_stream)) {
+      if (MAD_RECOVERABLE(mad_stream.error)) {
+        continue;
+      } else if (mad_stream.error == MAD_ERROR_BUFLEN) {
+        continue;
+      } else {
+        break;
+      }
+    }
+  }
+
+  mad_synth_finish(&mad_synth);
+  mad_frame_finish(&mad_frame);
+  mad_stream_finish(&mad_stream);
+
+  return output;
+}
+
+long MPEGFile::readFrame(const struct mad_pcm * pcm,
+               std::vector< std::vector<double> > & output,
+               long sampleOffset,
+               long startSample,
+               long endSample) {
+
+  // if the last sample in the frame is after our start sample
+  // and the first sample in the frame is before our end sample
+  if (((sampleOffset + pcm -> length) > startSample) and
+      (sampleOffset <= endSample)) {
+
+    for (short sampleNum = 0; sampleNum < pcm -> length; sampleNum ++) {
+      for (short channel = 0; channel < pcm -> channels; channel++) {
+
+        if ((sampleNum + sampleOffset >= startSample) and
+            (sampleNum + sampleOffset < endSample)) {
+
+          output[channel][sampleNum + sampleOffset - startSample] =
+            pcm -> samples[channel][sampleNum];
+
+        }
+      }
+    }
+  }
+
+  return sampleOffset + pcm -> length;
+}
diff --git a/Trash/MPEGFile.hpp b/Trash/MPEGFile.hpp
new file mode 100644
index 0000000..3e49a2a
--- /dev/null
+++ b/Trash/MPEGFile.hpp
@@ -0,0 +1,26 @@
+#ifndef MPEG_FILE_H
+#define MPEG_FILE_H
+
+#include <vector>
+
+class MPEGFile {
+  public:
+    /**
+     * Returns a vector whose size
+     * is the number of channels. Each channel
+     * contains the uncompressed audio data
+     * in the range 0...1.
+     */
+    static std::vector<std::vector<double> > read(
+        const std::string filePath,
+        double startSeconds,
+        double endSeoncds);
+
+    static long readFrame(const struct mad_pcm * pcm,
+                          std::vector< std::vector<double> > & ouput,
+                          long sampleOffset,
+                          long startSample,
+                          long endSample);
+};
+
+#endif
diff --git a/Trash/Tempo.cpp b/Trash/Tempo.cpp
new file mode 100644
index 0000000..b176afd
--- /dev/null
+++ b/Trash/Tempo.cpp
@@ -0,0 +1,267 @@
+#include <vector>
+#include <iostream>
+
+#include "Plotting/Plotting.hpp"
+#include "SampleSorter/Tempo.hpp"
+#include "SampleSorter/SpectralProcessing.hpp"
+
+double Tempo::tempoToSeconds(double tempo) {
+  return 1/tempo;
+}
+
+double Tempo::tempoToSamples(double tempo, long sampleRate) {
+  return secondsToSamples(tempoToSeconds(tempo), sampleRate);
+}
+
+double Tempo::secondsToSamples(double seconds, long sampleRate) {
+  return seconds * sampleRate;
+}
+
+double Tempo::samplesToSeconds(double samples, long sampleRate) {
+  return samples/double(sampleRate);
+}
+
+double Tempo::tempoToBins(double tempo, long hopSize, long sampleRate) {
+  return sampleRate/(tempo * hopSize);
+}
+
+double Tempo::binsToSamples(long bin, long hopSize) {
+  return bin * hopSize;
+}
+
+double Tempo::binsToSeconds(long bin, long hopSize, long sampleRate) {
+  return samplesToSeconds(binsToSamples(bin, hopSize), sampleRate);
+}
+
+
+double Tempo::tempoValue(double tempo, 
+                         double theOne,
+                         std::vector<double> onsets, 
+                         long hopSize,
+                         long sampleRate,
+                         bool bidirectional) {
+  double value = 0;
+  for (long i = 0; i < onsets.size(); i++) {
+    double distanceFromBeat = (i - theOne)/tempoToBins(tempo, hopSize, sampleRate);
+    if (bidirectional) {
+      distanceFromBeat = std::abs(distanceFromBeat - std::round(distanceFromBeat));
+    } else {
+      distanceFromBeat = distanceFromBeat - std::floor(distanceFromBeat);
+    }
+
+    value += distanceFromBeat * onsets[i];
+  }
+
+  return value;//tempoToBins(tempo, hopSize, sampleRate);
+}
+
+double Tempo::tempoValue(double tempo, 
+                         double theOne,
+                         std::vector<double> onsets, 
+                         long hopSize,
+                         long sampleRate,
+                         bool bidirectional) {
+  double value = 0;
+  for (long i = 0; i < onsets.size(); i++) {
+    double distanceFromBeat = (i - theOne)/tempoToBins(tempo, hopSize, sampleRate);
+    if (bidirectional) {
+      distanceFromBeat = std::abs(distanceFromBeat - std::round(distanceFromBeat));
+    } else {
+      distanceFromBeat = distanceFromBeat - std::floor(distanceFromBeat);
+    }
+
+    value += distanceFromBeat * onsets[i];
+  }
+
+  return value;//tempoToBins(tempo, hopSize, sampleRate);
+}
+
+
+double Tempo::correlationTempo(std::vector<double> onsets, 
+                               long hopSize, 
+                               long sampleRate) {
+  //Plotting::plotVector(onsets, hopSize);
+
+  // filter out frequencies less than 2 per clip
+  // Clips with tempo must contain at least 2 beats
+  double highPass = 2./binsToSeconds(onsets.size(), hopSize, sampleRate);
+  // Also beat cannot be below 1/8 persecond
+  highPass = std::max(highPass, 1.);
+
+  // No tempos will exist above 1000bpm
+  double lowPass = 1000/60.;
+
+  std::vector<double> correlation = SpectralProcessing::autoCorrelation(onsets,
+                                                    sampleRate/hopSize,
+                                                    highPass,
+                                                    lowPass);
+
+  // only take first half of correlation
+  long fftSize = correlation.size()/4 + 1;
+  fftw_complex * fft;
+  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
+  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(correlation.size()/2,
+                                              correlation.data(),
+                                              fft,
+                                              FFTW_ESTIMATE);
+
+  //Window the correlation
+  for (long i = 0; i < correlation.size()/2; i++) {
+    correlation[i] = correlation[i] * SpectralProcessing::hammingWindow(i, correlation.size()/2);
+  }
+
+  fftw_execute(fftPlan);
+  fftw_destroy_plan(fftPlan);
+
+  std::vector<std::pair<double, double> > peaks = 
+    SpectralProcessing::findPeaks(fft, fftSize, sampleRate/hopSize);
+
+  //Plotting::plotPair(peaks);
+
+  fftw_free(fft);
+
+  double max = 0;
+  double tempo = 0;
+  for (long i = 0; i < peaks.size(); i++) {
+    if (peaks[i].second > max) {
+      max = peaks[i].second;
+      tempo = peaks[i].first;
+    }
+  }
+
+  return tempo;
+}
+
+std::pair<double, double> Tempo::fineTuneTempo(double tempo, 
+    std::vector<double> onsets, 
+    double BPSrange, 
+    double BPSstepSize,
+    long hopSize,
+    long sampleRate) {
+
+  std::vector<double> tempoValues;
+
+  //Plotting::plotVector(onsets);
+
+  // first find the bidirectional minima
+  double minTempoValue = onsets.size();
+  double minTempo = tempo;
+  for (double tempoAdjustment = -BPSrange; 
+       tempoAdjustment <=BPSrange; 
+       tempoAdjustment += BPSstepSize) {
+    double trialTempo = tempo + tempoAdjustment;
+    // adjust it by size of beats to bins -> less bins -> easier
+    double thisMinValue = onsets.size();
+    double thisMinTempo = 0;
+    for (double theOne = 0; 
+         theOne < tempoToBins(trialTempo, hopSize, sampleRate); 
+         theOne ++) {
+      double newTempoValue = tempoValue(trialTempo, theOne, onsets, hopSize, sampleRate, true);
+      if (newTempoValue < minTempoValue) {
+        minTempoValue = newTempoValue;
+        minTempo = trialTempo;
+      }
+
+      if (newTempoValue < thisMinValue) {
+        thisMinValue = newTempoValue;
+        thisMinTempo = trialTempo;
+      }
+    }
+    tempoValues.push_back(thisMinValue);
+  }
+
+  // then search for the start of the beat
+  minTempoValue = onsets.size();
+  std::vector<double> minValues;
+  double minOne = 0;
+  for (double theOne = 0; 
+       theOne < tempoToBins(minTempo, hopSize, sampleRate)+1; 
+       theOne ++) {
+    double newTempoValue = tempoValue(minTempo, theOne, onsets, hopSize, sampleRate, false);
+    if (newTempoValue < minTempoValue) {
+      minTempoValue = newTempoValue;
+      minOne = theOne;
+    }
+    minValues.push_back(newTempoValue);
+  }
+
+  return std::make_pair(minTempo, binsToSeconds(minOne, hopSize, sampleRate));
+}
+
+std::vector<double> Tempo::fineTuneTempoFunction(
+    double tempo, 
+    size_t degrees,
+    std::vector<double> onsets, 
+    double errorPercentage,
+    double stepPercentage,
+    long hopSize,
+    long sampleRate) {
+
+  std::vector<double> function(degrees);
+  double minTempoValue = onsets.size();
+
+  // first find 0 degree. 
+  for (double tempoAdjustment = BPSPercentage * tempo; 
+       tempoAdjustment <=BPSrange; 
+       tempoAdjustment += BPSstepSize) {
+    double trialTempo = tempo + tempoAdjustment;
+    for (double theOne = 0; 
+         theOne < tempoToBins(trialTempo, hopSize, sampleRate); 
+         theOne ++) {
+      double newTempoValue = tempoValue(trialTempo, theOne, onsets, hopSize, sampleRate, true);
+      if (newTempoValue < minTempoValue) {
+        minTempoValue = newTempoValue;
+        minTempo = trialTempo;
+      }
+    }
+  }
+
+  for (size_t degree = 1; degree < degrees; degree ++) {
+    // for all possible values
+    // for all possible values of the one
+    // find the newTempoValue
+  }
+
+  // first find the bidirectional minima
+  double minTempoValue = onsets.size();
+  double minTempo = tempo;
+  for (double tempoAdjustment = BPSPercentage * tempo; 
+       tempoAdjustment <=BPSrange; 
+       tempoAdjustment += BPSstepSize) {
+    // adjust it by size of beats to bins -> less bins -> easier
+    double thisMinValue = onsets.size();
+    double thisMinTempo = 0;
+    for (double theOne = 0; 
+         theOne < tempoToBins(trialTempo, hopSize, sampleRate); 
+         theOne ++) {
+      double newTempoValue = tempoValue(trialTempo, theOne, onsets, hopSize, sampleRate, true);
+      if (newTempoValue < minTempoValue) {
+        minTempoValue = newTempoValue;
+        minTempo = trialTempo;
+      }
+
+      if (newTempoValue < thisMinValue) {
+        thisMinValue = newTempoValue;
+        thisMinTempo = trialTempo;
+      }
+    }
+    tempoValues.push_back(thisMinValue);
+  }
+
+  // then search for the start of the beat
+  minTempoValue = onsets.size();
+  std::vector<double> minValues;
+  double minOne = 0;
+  for (double theOne = 0; 
+       theOne < tempoToBins(minTempo, hopSize, sampleRate)+1; 
+       theOne ++) {
+    double newTempoValue = tempoValue(minTempo, theOne, onsets, hopSize, sampleRate, false);
+    if (newTempoValue < minTempoValue) {
+      minTempoValue = newTempoValue;
+      minOne = theOne;
+    }
+    minValues.push_back(newTempoValue);
+  }
+
+  return std::make_pair(minTempo, binsToSeconds(minOne, hopSize, sampleRate));
+}
diff --git a/Trash/Tempo.hpp b/Trash/Tempo.hpp
new file mode 100644
index 0000000..c5bafc6
--- /dev/null
+++ b/Trash/Tempo.hpp
@@ -0,0 +1,35 @@
+#ifndef TEMPO_H
+#define TEMPO_H
+
+class Tempo {
+  public:
+
+  static double tempoToSeconds(double tempo);
+  static double tempoToSamples(double tempo, long sampleRate);
+  static double secondsToSamples(double seconds, long sampleRate);
+  static double samplesToSeconds(double samples, long sampleRate);
+  static double tempoToBins(double tempo, long hopSize, long sampleRate);
+  static double binsToSamples(long bin, long hopSize);
+  static double binsToSeconds(long bin, long hopSize, long sampleRate);
+
+  static double correlationTempo(std::vector<double> offsets,
+                                 long hopSize,
+                                 long sampleRate);
+
+  static double tempoValue(double tempo,
+                           double theOne,
+                           std::vector<double> offsets,
+                           long hopSize,
+                           long sampleRate,
+                           bool bidirectional);
+
+
+  static std::pair<double, double> fineTuneTempo(double tempo,
+                              std::vector<double> offsets,
+                              double BPSrange,
+                              double BPSstepSize,
+                              long hopSize,
+                              long sampleRate);
+};
+
+#endif
diff --git a/Trash/TempoFunction.cpp b/Trash/TempoFunction.cpp
new file mode 100644
index 0000000..ddbdb6a
--- /dev/null
+++ b/Trash/TempoFunction.cpp
@@ -0,0 +1,408 @@
+#include <vector>
+#include <cmath>
+#include <random>
+
+#include <fftw3.h>
+
+#include "SampleSorter/TempoFunction.hpp"
+#include "SampleSorter/Units.hpp"
+#include "SampleSorter/SpectralProcessing.hpp"
+#include "Plotting/Plotting.hpp"
+
+TempoFunction::TempoFunction() {
+  avgTempo = 1;
+  totalSeconds = 0;
+}
+
+double TempoFunction::getTotalSeconds() const {
+  return totalSeconds;
+}
+
+double TempoFunction::getTheOne() const {
+  return aCoefficients[0]/2.;
+}
+
+TempoFunction::TempoFunction(double avgTempo_,
+              double theOne,
+              double totalSeconds_) {
+  aCoefficients.resize(1);
+  bCoefficients.resize(1);
+  avgTempo = avgTempo_;
+  aCoefficients[0] = theOne * 2.;
+  totalSeconds = totalSeconds_;
+}
+
+void TempoFunction::fineTuneAvgTempo(
+    double guessTempo,
+    double percentageError,
+    double percentageStep,
+    const std::vector<double> & onsets,
+    long hopSize,
+    long sampleRate
+    ) {
+
+  double minValue = onsets.size();
+  double minTempo = guessTempo;
+  double minOne = 0;
+
+  for (double trialTempo = guessTempo/2.;
+       trialTempo < guessTempo * 2.;
+       trialTempo *= std::pow(2, 0.001)) {
+
+    avgTempo = trialTempo;
+    aCoefficients[0] = 0;
+
+    double theFirstBeat = getKthBeat(1);
+
+    for (double theOne = 0;
+        theOne < Units::secondsToBins(theFirstBeat, hopSize, sampleRate);
+        theOne ++) {
+
+      aCoefficients[0] = theOne;
+
+      double newValue = getValue(onsets, hopSize, sampleRate);
+      if (newValue < minValue) {
+        minValue = newValue;
+        minTempo = trialTempo;
+        minOne = theOne;
+      }
+    }
+  }
+
+  avgTempo = minTempo;
+  aCoefficients[0] = minOne;
+}
+
+void TempoFunction::fineTuneTheOne(
+    const std::vector<double> & onsets,
+    long hopSize,
+    long sampleRate) {
+  // set the one to zero
+  aCoefficients[0] = 0;
+
+  // get the first beat
+  double beatOne = getKthBeat(1);
+
+  double minValue = onsets.size(), minOne = 0;
+
+  // for all possible ones
+  for (double onePosition = 0;
+      onePosition < beatOne;
+      onePosition += 0.001) {
+    aCoefficients[0] = onePosition;
+    double value = 0;
+    for (long bin = 0; bin < onsets.size(); bin ++) {
+      double binSeconds = Units::binsToSeconds(bin, hopSize, sampleRate);
+      value += onsets[bin] * std::fmod(getBeatNum(binSeconds), 1.);
+    }
+    if (value < minValue) {
+      minValue = value;
+      minOne = onePosition;
+    }
+  }
+  aCoefficients[0] = minOne;
+}
+
+double TempoFunction::getAvgTempo() const {
+  return avgTempo;
+}
+
+size_t TempoFunction::getDegree() const {
+  return aCoefficients.size();
+}
+
+double TempoFunction::getBeatNum(double time) const {
+  // actual tempo
+  double beatNum = avgTempo * time;
+  // constant offset
+  beatNum += aCoefficients[0]/2.;
+
+  // Fourier series approximation
+  for (size_t i = 1; i < getDegree(); i++) {
+    beatNum += aCoefficients[i] * 
+               std::cos((2.*M_PI*i*time)/totalSeconds);
+    beatNum += bCoefficients[i] * 
+               std::sin((2.*M_PI*i*time)/totalSeconds);
+  }
+
+  return beatNum;
+}
+
+// The derivative
+double TempoFunction::getTempo(double time) const {
+  double tempo = avgTempo;
+
+  // Fourier series approximation
+  for (size_t i = 1; i < getDegree(); i++) {
+    tempo -= aCoefficients[i] * 
+             ((2*M_PI*i)/totalSeconds) *
+             std::sin(2.*M_PI*i*time/totalSeconds);
+
+    tempo += bCoefficients[i] * 
+             ((2*M_PI*i)/totalSeconds) *
+             std::cos(2.*M_PI*i*time/totalSeconds);
+  }
+
+  return tempo;
+}
+
+double TempoFunction::distanceFromBeat(double time) const {
+  return 1 - std::cos(2.*M_PI*getBeatNum(time));
+}
+
+double TempoFunction::getValue(
+    const std::vector<double> & onsets,
+    long hopSize,
+    long sampleRate
+    ) const {
+  double value = 0;
+  for (size_t bin = 0; bin < onsets.size(); bin ++) {
+    double time = Units::binsToSeconds(bin, hopSize, sampleRate);
+    value += distanceFromBeat(time) * onsets[bin];
+  }
+  
+  return value;
+}
+
+void TempoFunction::getBinGradient(
+    size_t bin,
+    std::vector<double> & aNabla,
+    std::vector<double> & bNabla,
+    const std::vector<double> & onsets,
+    long hopSize,
+    long sampleRate
+    ) const {
+
+  double time = Units::binsToSeconds(bin,hopSize, sampleRate);
+
+  double multiplier = onsets[bin];
+  multiplier *= 2.*M_PI * std::sin(2.*M_PI*getBeatNum(time));
+
+  aNabla[0] = 0.5 * multiplier;
+
+  for (size_t i = 1; i < getDegree(); i++) {
+    aNabla[i] = std::cos(2.*M_PI*i*time/totalSeconds);
+    aNabla[i] *= multiplier;
+
+    bNabla[i] = std::sin(2.*M_PI*i*time/totalSeconds);
+    bNabla[i] *= multiplier;
+  }
+}
+
+void TempoFunction::gradientDescent(
+    const std::vector<double> & onsets,
+    long hopSize, 
+    long sampleRate
+    ) {
+
+  std::vector<double> aNabla(getDegree());
+  std::vector<double> bNabla(getDegree());
+
+  std::vector<double> aDelta(getDegree()); 
+  std::fill(aDelta.begin(), aDelta.end(), 0);
+  std::vector<double> bDelta(getDegree()); 
+  std::fill(bDelta.begin(), bDelta.end(), 0);
+
+  double learningRate;
+  double mass = 0.5;
+
+  std::srand(std::time(0));
+
+  long numIterations = onsets.size() * 4;
+  //std::cout << "numIterations" << numIterations << std::endl;
+  std::vector<double> values(numIterations);
+
+  for (long trial = 0; trial < numIterations; trial++) {
+    learningRate = 0.000001*(numIterations - trial)/double(numIterations);
+
+    // choose a random bin
+    size_t bin = std::rand() % onsets.size();
+
+    getBinGradient(bin, aNabla, bNabla, onsets, hopSize, sampleRate);
+
+    for (size_t i = 0; i < getDegree(); i++) {
+      aDelta[i] = learningRate * aNabla[i] + mass * aDelta[i];
+      aCoefficients[i] = aCoefficients[i] - aDelta[i];
+
+      bDelta[i] = learningRate * bNabla[i] + mass * bDelta[i];
+      bCoefficients[i] = bCoefficients[i] - bDelta[i];
+    }
+
+    values[trial] = getValue(onsets, hopSize, sampleRate);
+  }
+
+}
+
+
+double TempoFunction::getKthBeat(double k, double accuracy) const {
+  //std::cout << "starting" << k << std::endl;
+  double guess = k/avgTempo + aCoefficients[0];
+
+  double error = getBeatNum(guess) - k;
+  while (std::abs(error) > accuracy) {
+    guess -= error/getTempo(guess);
+    //std::cout << "tempo: " << getTempo(guess) << std::endl;
+    //std::cout << "error: " << error << std::endl;
+    //std::cout << "guess: " << guess << std::endl;
+    error = getBeatNum(guess) - k;
+  }
+
+  //std::cout << "ending" << k << std::endl;
+
+  return guess;
+}
+
+TempoFunction::TempoFunction(
+    std::vector< std::vector<double> > & audio,
+    long sampleRate,
+    size_t degrees, 
+    double percentageError,
+    double percentageStep) {
+
+  std::cout << audio.size() << ", " << audio[0].size() << std::endl;
+  std::cout << sampleRate << std::endl;
+
+  //std::cout << "starting tempo" << std::endl;
+
+  totalSeconds = Units::samplesToSeconds(audio[0].size(), sampleRate);
+
+  long hopSize = 1024;
+  long windowRatio = 2;
+  std::vector<double> onsets =
+    SpectralProcessing::onsetEnergy(audio, hopSize, windowRatio);
+
+  double guessTempo = correlationTempo(onsets, hopSize, sampleRate);
+
+  std::cout << "guessed tempo: " << guessTempo * 60. << std::endl;
+
+  aCoefficients.resize(1);
+  bCoefficients.resize(1);
+  std::fill(aCoefficients.begin(), aCoefficients.end(), 0);
+  std::fill(bCoefficients.begin(), bCoefficients.end(), 0);
+
+  fineTuneAvgTempo(guessTempo, 
+                   percentageError, 
+                   percentageStep, 
+                   onsets, 
+                   hopSize, 
+                   sampleRate);
+
+  std::cout << "fine tuned tempo: " << getAvgTempo() * 60. << std::endl;
+
+  //aCoefficients.resize(degrees);
+  //bCoefficients.resize(degrees);
+  //std::fill(aCoefficients.begin(), aCoefficients.end(), 0);
+  //std::fill(bCoefficients.begin(), bCoefficients.end(), 0);
+  //gradientDescent(onsets, hopSize, sampleRate);
+
+  fineTuneTheOne(onsets, hopSize, sampleRate);
+
+  //for (long i = 0; i < getDegree(); i++) {
+    //std::cout << "a[" << i << "]=" << aCoefficients[i] << std::endl;
+    //std::cout << "a[" << i << "]=" << bCoefficients[i] << std::endl;
+  //}
+
+  //plotOnsetsWithBeats(onsets, hopSize, sampleRate);
+  //plotBeats(0.001);
+  //plotTempo(0.001);
+  //
+
+}
+
+double TempoFunction::correlationTempo(
+                     const std::vector<double> & onsets, 
+                     long hopSize, 
+                     long sampleRate) {
+
+  // filter out frequencies less than 2 per clip
+  // Clips with tempo must contain at least 2 beats
+  double highPass = 2./Units::binsToSeconds(onsets.size(), hopSize, sampleRate);
+  // Also beat cannot be below 1/8 persecond
+  highPass = std::max(highPass, 1.);
+
+  // No tempos will exist above 1000bpm
+  double lowPass = 1000/60.;
+
+  std::vector<double> correlation = SpectralProcessing::autoCorrelation(onsets,
+                                                    sampleRate/hopSize,
+                                                    highPass,
+                                                    lowPass);
+
+  // only take first half of correlation
+  long fftSize = correlation.size()/4 + 1;
+  fftw_complex * fft;
+  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
+  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(correlation.size()/2,
+                                              correlation.data(),
+                                              fft,
+                                              FFTW_ESTIMATE);
+
+  //Window the correlation
+  for (long i = 0; i < correlation.size()/2; i++) {
+    correlation[i] = correlation[i] 
+              * SpectralProcessing::hammingWindow(i, correlation.size()/2);
+  }
+
+  fftw_execute(fftPlan);
+  fftw_destroy_plan(fftPlan);
+
+  std::vector<std::pair<double, double> > peaks = 
+    SpectralProcessing::findPeaks(fft, fftSize, sampleRate/hopSize);
+
+  //Plotting::plotPair(peaks);
+
+  fftw_free(fft);
+
+  double max = 0;
+  double tempo = 0;
+  for (long i = 0; i < peaks.size(); i++) {
+    if (peaks[i].second > max) {
+      max = peaks[i].second;
+      tempo = peaks[i].first;
+    }
+  }
+
+  return tempo;
+}
+
+void TempoFunction::plotTempo(double resolution) const {
+  std::vector<double> tempo(totalSeconds/resolution);
+  for (long i = 0; i*resolution < totalSeconds; i++) {
+    tempo[i] = 60.*getTempo(i*resolution);
+  }
+  Plotting::plotVector(tempo);
+}
+
+void TempoFunction::plotBeats(double resolution) const {
+  std::vector<double> beats(totalSeconds/resolution);
+  for (long i = 0; i*resolution < totalSeconds; i++) {
+    beats[i] = getBeatNum(i*resolution);
+  }
+  Plotting::plotVector(beats);
+}
+
+void TempoFunction::plotOnsetsWithBeats(
+    const std::vector<double> & onsets,
+    long hopSize,
+    long sampleRate
+    ) const {
+
+  std::vector<double> beats;
+  long k = 0;
+  double beat;
+  do {
+    beat = getKthBeat(k);
+    beats.push_back(beat);
+    k += 1;
+  } while (beat < totalSeconds);
+
+  std::vector<std::pair<double, double> > onsetSeconds(onsets.size());
+
+  for (long bin = 0; bin < onsets.size(); bin++) {
+    onsetSeconds[bin] = std::make_pair(Units::binsToSeconds(bin, hopSize, sampleRate), onsets[bin]);
+  }
+
+  Plotting::plotLineAndMarkers(onsetSeconds, beats, 0.5);
+}
+
+// plot onsets with beats
diff --git a/Trash/TempoFunction.hpp b/Trash/TempoFunction.hpp
new file mode 100644
index 0000000..16b235e
--- /dev/null
+++ b/Trash/TempoFunction.hpp
@@ -0,0 +1,92 @@
+#ifndef TEMPO_FUNCTION_H
+#define TEMPO_FUNCTION_H
+
+#include <vector>
+#include <cmath>
+#include <random>
+
+class TempoFunction {
+  private:
+    double avgTempo;
+    double totalSeconds;
+    std::vector<double> aCoefficients;
+    std::vector<double> bCoefficients;
+  public:
+    TempoFunction();
+    double getTotalSeconds() const;
+    double getTheOne() const;
+
+    TempoFunction(double avgTempo_,
+                  double theOne,
+                  double totalSeconds_);
+
+    TempoFunction(
+        std::vector<std::vector<double> > & audio,
+        long sampleRate,
+        size_t degrees, 
+        double percentageError,
+        double percentageStep);
+
+    void fineTuneAvgTempo(
+        double guessTempo,
+        double percentageError,
+        double percentageStep,
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate);
+
+    void gradientDescent(
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate
+        );
+
+    void fineTuneTheOne(
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate
+        );
+
+    double getAvgTempo() const;
+
+    size_t getDegree() const;
+
+    double getBeatNum(double time) const;
+
+    // The derivative
+    double getTempo(double time) const;
+
+    double distanceFromBeat(double time) const;
+
+    double getValue(
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate
+        ) const;
+
+    void getBinGradient(
+        size_t bin,
+        std::vector<double> & aNabla,
+        std::vector<double> & bNabla,
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate
+        ) const;
+
+    double getKthBeat(double k, double accuracy = 0.005) const;
+
+    static double correlationTempo(
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate);
+
+    void plotBeats(double resolution) const;
+    void plotTempo(double resolution) const;
+    void plotOnsetsWithBeats(
+        const std::vector<double> & onsets,
+        long hopSize,
+        long sampleRate) const;
+    
+};
+
+#endif
diff --git a/include/Plotting/Plotting.hpp b/include/Plotting/Plotting.hpp
index 7c81126..f137841 100644
--- a/include/Plotting/Plotting.hpp
+++ b/include/Plotting/Plotting.hpp
@@ -7,6 +7,10 @@ class Plotting {
   public:
     static void plotVector(std::vector<double> y, double xAxisScale = 1, double xAxisOffset = 0);
     static void plotPair(std::vector<std::pair<double, double> > xy);
+    static void plotLineAndMarkers(
+        std::vector<std::pair<double, double> > line,
+        std::vector<double> markers,
+        double pointHeight = 1);
 };
 
 #endif
diff --git a/include/SampleSorter/AbletonSampleFile.hpp b/include/SampleSorter/AbletonSampleFile.hpp
index 884dfb3..f08ee52 100644
--- a/include/SampleSorter/AbletonSampleFile.hpp
+++ b/include/SampleSorter/AbletonSampleFile.hpp
@@ -5,10 +5,17 @@
 
 #include <tinyxml2.h>
 
-#include <sndfile.hh>
-
 #include "SampleSorter/SampleFile.hpp"
 
+//class AbletonNode {
+  //private:
+  //public:
+    //// constructor with list
+    //AbletonNode
+    
+
+//}
+
 class AbletonSampleFile : public SampleFile {
   private:
     tinyxml2::XMLDocument doc;
@@ -17,13 +24,41 @@ class AbletonSampleFile : public SampleFile {
     double startSeconds;
     double endSeconds;
 
-    SndfileHandle audioFile;
-
+    // opens file
     void getDoc();
-    void readDoc();
 
-    std::vector< std::vector<double> > getWaves();
-    long getSampleRate() const;
+    // parses file
+    // returns true iff preprocessed
+    bool readDoc();
+
+    // writes file
+    void setDoc();
+
+    virtual std::vector< std::vector<double> > extractAudio(long * sampleRate);
+
+    // opens and parses file
+    virtual bool readMetaData();
+
+    // Constants
+    //static const std::string VALUE = "Value";
+
+    // Document reading functions
+    tinyxml2::XMLElement * getNode(std::initializer_list<std::string> nodes);
+    tinyxml2::XMLElement * getAudioNode();
+    tinyxml2::XMLElement * getLoopNode();
+    tinyxml2::XMLElement * getSortDataNode();
+
+    // Direct ones
+    tinyxml2::XMLElement * getLoopStartNode();
+    tinyxml2::XMLElement * getLoopEndNode();
+    tinyxml2::XMLElement * getHiddenLoopStartNode();
+    tinyxml2::XMLElement * getHiddenLoopEndNode();
+    tinyxml2::XMLElement * getNameNode();
+    tinyxml2::XMLElement * getPitchFineNode();
+
+    
+    // various fetchers
+    std::string fetchReferenceFilePath();
 
   public:
     AbletonSampleFile(std::string filePath);
@@ -33,7 +68,9 @@ class AbletonSampleFile : public SampleFile {
     std::string getReferenceFilePath() const;
     std::string getReferenceFileName() const;
 
-    double getSampleLength() const;
+    virtual double getSampleSeconds() const;
+
+    void writeToFile();
 };
 
 #endif
diff --git a/include/SampleSorter/AudioFile.hpp b/include/SampleSorter/AudioFile.hpp
new file mode 100644
index 0000000..efaa6b9
--- /dev/null
+++ b/include/SampleSorter/AudioFile.hpp
@@ -0,0 +1,26 @@
+#ifndef AUDIO_FILE_H
+#define AUDIO_FILE_H
+
+#include <vector>
+#include <string>
+
+#include <libavcodec/avcodec.h>
+
+class AudioFile {
+  public:
+    static std::vector<std::vector<double> > read(
+        std::string fileName,
+        double startSeconds,
+        double endSeconds,
+        long * sampleRate);
+
+    static void readFrame(
+        const AVCodecContext * codecContext,
+        const AVFrame * frame,
+        std::vector<std::vector<double> > & output,
+        long startSample,
+        long endSample,
+        long sampleOffset);
+};
+
+#endif
diff --git a/include/SampleSorter/AudioSample.hpp b/include/SampleSorter/AudioSample.hpp
index d897531..de8184d 100644
--- a/include/SampleSorter/AudioSample.hpp
+++ b/include/SampleSorter/AudioSample.hpp
@@ -4,21 +4,30 @@
 #include <vector>
 
 #include "SampleSorter/Octave.hpp"
+#include "SampleSorter/Tempo.hpp"
 
 class AudioSample {
   private:
     long tuningCents;
-    double tempo;
-    double theOne;
+    Tempo tempo;
     std::vector<Octave> chords;
+    double totalSeconds;
+    long sampleRate;
 
-    void tune(std::vector<std::vector<double> > audio, long sampleRate);
-    void findBeat(std::vector<std::vector<double> > audio, long sampleRate);
-    void findChords(std::vector<std::vector<double> > audio, long sampleRate);
+    void tune(std::vector<std::vector<double> > & audio);
+    void findBeat(std::vector<std::vector<double> > & audio);
+    void findChords(std::vector<std::vector<double> > & audio);
   public:
     AudioSample();
-    AudioSample(std::vector<std::vector<double> > audio, long sampleRate);
+    AudioSample(long tuningCents_,
+                double rawBeat,
+                double theOne,
+                double totalSeconds,
+                long sampleRate,
+                std::vector<Octave> chords);
+    AudioSample(std::vector<std::vector<double> > & audio, long _sampleRate);
 
+    double getTotalSeconds() const;
     long getTuningCents() const;
     double getTuningCentsFreqRatio() const;
     double getBeatRaw() const;
@@ -26,6 +35,7 @@ class AudioSample {
     double getTheOneRaw() const;
     double getTheOneWithTuning() const;
     std::vector<Octave> getChords() const;
+    long getSampleRate() const;
 };
 
 #endif
diff --git a/include/SampleSorter/ProcessingException.hpp b/include/SampleSorter/ProcessingException.hpp
new file mode 100644
index 0000000..6612a3d
--- /dev/null
+++ b/include/SampleSorter/ProcessingException.hpp
@@ -0,0 +1,15 @@
+#include <string>
+#include <exception>
+
+class ProcessingException: public std::exception {
+  private:
+    const char * errorMessage;
+  public:
+    ProcessingException(std::string errorMessage_) {
+      errorMessage = errorMessage_.c_str();
+    }
+
+    virtual const std::string getMessage() const throw() {
+      return errorMessage;
+    }
+};
diff --git a/include/SampleSorter/SampleFile.hpp b/include/SampleSorter/SampleFile.hpp
index ac015b3..521e958 100644
--- a/include/SampleSorter/SampleFile.hpp
+++ b/include/SampleSorter/SampleFile.hpp
@@ -7,30 +7,20 @@
 
 class SampleFile {
   private:
-    AudioSample sample;
-
-    virtual std::vector< std::vector<double> > getWaves() = 0;
-    virtual long getSampleRate() const = 0;
+    virtual std::vector< std::vector<double> > extractAudio(long * sampleRate) = 0;
+    virtual bool readMetaData() = 0;
   protected:
     std::string filePath;
+    AudioSample sample;
   public:
     SampleFile(std::string filePath_);
-    void process();
+    bool process();
 
     std::string getFilePath();
     std::string getFileName();
 
-    virtual double getSampleLength() const = 0;
+    virtual double getSampleSeconds() const = 0;
     AudioSample * getAudioSample();
 };
 
-//extern "C" {
-  //SampleFile * SampleFile(char * fileName) {
-    //return new SampleFile(fileName);
-  //}
-  //double getTempo(SampleFile * s) {
-    //return s -> getAudioSample -> getTempo();
-  //}
-//}
-
 #endif
diff --git a/include/SampleSorter/Tempo.hpp b/include/SampleSorter/Tempo.hpp
index c5bafc6..649ccc2 100644
--- a/include/SampleSorter/Tempo.hpp
+++ b/include/SampleSorter/Tempo.hpp
@@ -1,35 +1,106 @@
-#ifndef TEMPO_H
-#define TEMPO_H
+#ifndef TEMPO_FUNCTION_H
+#define TEMPO_FUNCTION_H
+
+#include <vector>
+#include <cmath>
+#include <random>
 
 class Tempo {
+  private:
+    /**
+     * The average tempo of the
+     * sample
+     */
+    double tempo;
+    /**
+     * The starting beat of the
+     * sample measured in bins
+     */
+    double theOneBin;
+
+    static const long HOP_SIZE;
+    long sampleRate;
   public:
+    /**
+     * The tempo.
+     */
+    double getTempo() const;
+    /**
+     * Returns the starting
+     * beat in bins
+     */
+    double getTheOneBin() const;
+    /**
+     * Returns the starting beat
+     * in seconds
+     */
+    double getTheOne() const;
+
+    Tempo();
+    Tempo(double tempo_, double theOneBin_, long sampleRate);
+
+    /**
+     * Takes an audio file
+     * and measures the
+     * tempo and the one
+     */
+    Tempo(
+        const std::vector<std::vector<double> > & audio,
+        long sampleRate,
+        double percentageError,
+        int tempoSteps,
+        int oneSteps
+        );
+
+    /**
+     * For a given guess tempo,
+     * this finds the tempo whose length
+     * in seconds is within the percentage
+     * error bounds that maximizes the energy
+     * around beats. This energy is measured
+     * by the getValue function.
+     */
+    void fineTuneTempo(
+        const double percentageError,
+        const int steps,
+        const std::vector<double> & onsets
+        );
+
+    /**
+     * Fine tunes the one so that most of the energy
+     * lies after each beat.
+     */
+    void fineTuneTheOne(
+        const std::vector<double> & onsets,
+        const int steps
+        );
+
+
+    /**
+     * Computes the distance from a beat.
+     * If bidirectional, the distance is,
+     * from either the previous or next beat,
+     * whichever is closer.
+     * If not bidirectional, the distance is
+     * the distance from the previous beat.
+     */
+    double distanceFromBeat(double bin, bool bidirectional, long hopSize, long sampleRate) const;
+
+    double getValue(
+        const std::vector<double> & onsets,
+        bool bidirectional
+        ) const;
+
+    void findCorrelationTempo(
+        const std::vector<double> & onsets
+        );
 
-  static double tempoToSeconds(double tempo);
-  static double tempoToSamples(double tempo, long sampleRate);
-  static double secondsToSamples(double seconds, long sampleRate);
-  static double samplesToSeconds(double samples, long sampleRate);
-  static double tempoToBins(double tempo, long hopSize, long sampleRate);
-  static double binsToSamples(long bin, long hopSize);
-  static double binsToSeconds(long bin, long hopSize, long sampleRate);
-
-  static double correlationTempo(std::vector<double> offsets,
-                                 long hopSize,
-                                 long sampleRate);
-
-  static double tempoValue(double tempo,
-                           double theOne,
-                           std::vector<double> offsets,
-                           long hopSize,
-                           long sampleRate,
-                           bool bidirectional);
-
-
-  static std::pair<double, double> fineTuneTempo(double tempo,
-                              std::vector<double> offsets,
-                              double BPSrange,
-                              double BPSstepSize,
-                              long hopSize,
-                              long sampleRate);
+    void plotBeats(double resolution) const;
+    void plotTempo(double resolution) const;
+    void plotOnsetsWithBeats(
+        const std::vector<double> & onsets
+        ) const;
+    
 };
 
 #endif
diff --git a/include/SampleSorter/Units.hpp b/include/SampleSorter/Units.hpp
new file mode 100644
index 0000000..1b6bad5
--- /dev/null
+++ b/include/SampleSorter/Units.hpp
@@ -0,0 +1,33 @@
+#ifndef UNITS_H
+#define UNITS_H
+
+#include <vector>
+#include <iostream>
+
+class Units {
+  public:
+
+    static double tempoToSeconds(double tempo);
+
+    static double tempoToSamples(double tempo, long sampleRate);
+
+    static double secondsToSamples(double seconds, long sampleRate);
+
+    static double secondsToTempo(double seconds);
+
+    static double samplesToSeconds(double samples, long sampleRate);
+
+    static double tempoToBins(double tempo, long hopSize, long sampleRate);
+
+    static double binsToSamples(double bin, long hopSize);
+
+    static double binsToSeconds(double bin, long hopSize, long sampleRate);
+
+    static double binsToTempo(double bin, long hopSize, long sampleRate);
+
+    static double samplesToBins(double samples, long hopSize);
+
+    static double secondsToBins(double seconds, long hopSize, long sampleRate);
+};
+
+#endif
diff --git a/modules/FindFFmpeg.cmake b/modules/FindFFmpeg.cmake
new file mode 100644
index 0000000..bd42350
--- /dev/null
+++ b/modules/FindFFmpeg.cmake
@@ -0,0 +1,147 @@
+# Locate ffmpeg
+# This module defines
+# FFMPEG_LIBRARIES
+# FFMPEG_FOUND, if false, do not try to link to ffmpeg
+# FFMPEG_INCLUDE_DIR, where to find the headers
+#
+# $FFMPEG_DIR is an environment variable that would
+# correspond to the ./configure --prefix=$FFMPEG_DIR
+#
+# Created by Robert Osfield.
+
+
+#In ffmpeg code, old version use "#include <header.h>" and newer use "#include <libname/header.h>"
+#In OSG ffmpeg plugin, we use "#include <header.h>" for compatibility with old version of ffmpeg
+
+#We have to search the path which contain the header.h (usefull for old version)
+#and search the path which contain the libname/header.h (usefull for new version)
+
+#Then we need to include ${FFMPEG_libname_INCLUDE_DIRS} (in old version case, use by ffmpeg header and osg plugin code)
+#                                                       (in new version case, use by ffmpeg header) 
+#and ${FFMPEG_libname_INCLUDE_DIRS/libname}             (in new version case, use by osg plugin code)
+
+
+# Macro to find header and lib directories
+# example: FFMPEG_FIND(AVFORMAT avformat avformat.h)
+MACRO(FFMPEG_FIND varname shortname headername)
+    # old version of ffmpeg put header in $prefix/include/[ffmpeg]
+    # so try to find header in include directory
+    FIND_PATH(FFMPEG_${varname}_INCLUDE_DIRS ${headername}
+        PATHS
+        ${FFMPEG_ROOT}/include/lib${shortname}
+        $ENV{FFMPEG_DIR}/include/lib${shortname}
+        ~/Library/Frameworks/lib${shortname}
+        /Library/Frameworks/lib${shortname}
+        /usr/local/include/lib${shortname}
+        /usr/include/lib${shortname}
+        /sw/include/lib${shortname} # Fink
+        /opt/local/include/lib${shortname} # DarwinPorts
+        /opt/csw/include/lib${shortname} # Blastwave
+        /opt/include/lib${shortname}
+        /usr/freeware/include/lib${shortname}
+        PATH_SUFFIXES ffmpeg
+        DOC "Location of FFMPEG Headers"
+    )
+
+    FIND_PATH(FFMPEG_${varname}_INCLUDE_DIRS ${headername}
+        PATHS
+        ${FFMPEG_ROOT}/include
+        $ENV{FFMPEG_DIR}/include
+        ~/Library/Frameworks
+        /Library/Frameworks
+        /usr/local/include
+        /usr/include
+        /sw/include # Fink
+        /opt/local/include # DarwinPorts
+        /opt/csw/include # Blastwave
+        /opt/include
+        /usr/freeware/include
+        PATH_SUFFIXES ffmpeg
+        DOC "Location of FFMPEG Headers"
+    )
+
+    FIND_LIBRARY(FFMPEG_${varname}_LIBRARIES
+        NAMES ${shortname}
+        PATHS
+        ${FFMPEG_ROOT}/lib
+        $ENV{FFMPEG_DIR}/lib
+        ~/Library/Frameworks
+        /Library/Frameworks
+        /usr/local/lib
+        /usr/local/lib64
+        /usr/lib
+        /usr/lib64
+        /sw/lib
+        /opt/local/lib
+        /opt/csw/lib
+        /opt/lib
+        /usr/freeware/lib64
+        DOC "Location of FFMPEG Libraries"
+    )
+
+    IF (FFMPEG_${varname}_LIBRARIES AND FFMPEG_${varname}_INCLUDE_DIRS)
+        SET(FFMPEG_${varname}_FOUND 1)
+    ENDIF(FFMPEG_${varname}_LIBRARIES AND FFMPEG_${varname}_INCLUDE_DIRS)
+
+ENDMACRO(FFMPEG_FIND)
+
+SET(FFMPEG_ROOT "$ENV{FFMPEG_DIR}" CACHE PATH "Location of FFMPEG")
+
+# find stdint.h
+IF(WIN32)
+
+    FIND_PATH(FFMPEG_STDINT_INCLUDE_DIR stdint.h
+        PATHS
+        ${FFMPEG_ROOT}/include
+        $ENV{FFMPEG_DIR}/include
+        ~/Library/Frameworks
+        /Library/Frameworks
+        /usr/local/include
+        /usr/include
+        /sw/include # Fink
+        /opt/local/include # DarwinPorts
+        /opt/csw/include # Blastwave
+        /opt/include
+        /usr/freeware/include
+        PATH_SUFFIXES ffmpeg
+        DOC "Location of FFMPEG stdint.h Header"
+    )
+
+    IF (FFMPEG_STDINT_INCLUDE_DIR)
+        SET(STDINT_OK TRUE)
+    ENDIF()
+
+ELSE()
+
+    # SET(STDINT_OK TRUE)
+
+ENDIF()
+
+FFMPEG_FIND(LIBAVFORMAT avformat avformat.h)
+FFMPEG_FIND(LIBAVDEVICE avdevice avdevice.h)
+FFMPEG_FIND(LIBAVCODEC  avcodec  avcodec.h)
+FFMPEG_FIND(LIBAVUTIL   avutil   avutil.h)
+FFMPEG_FIND(LIBSWSCALE  swscale  swscale.h)  # not sure about the header to look for here.
+
+SET(FFMPEG_FOUND "NO")
+# Note we don't check FFMPEG_LIBSWSCALE_FOUND here, it's optional.
+IF   (FFMPEG_LIBAVFORMAT_FOUND AND FFMPEG_LIBAVDEVICE_FOUND AND FFMPEG_LIBAVCODEC_FOUND AND FFMPEG_LIBAVUTIL_FOUND AND STDINT_OK)
+
+    SET(FFMPEG_FOUND "YES")
+
+    SET(FFMPEG_INCLUDE_DIRS ${FFMPEG_LIBAVFORMAT_INCLUDE_DIRS})
+
+    SET(FFMPEG_LIBRARY_DIRS ${FFMPEG_LIBAVFORMAT_LIBRARY_DIRS})
+
+    # Note we don't add FFMPEG_LIBSWSCALE_LIBRARIES here, it will be added if found later.
+    SET(FFMPEG_LIBRARIES
+        ${FFMPEG_LIBAVFORMAT_LIBRARIES}
+        ${FFMPEG_LIBAVDEVICE_LIBRARIES}
+        ${FFMPEG_LIBAVCODEC_LIBRARIES}
+        ${FFMPEG_LIBAVUTIL_LIBRARIES})
+
+ELSE ()
+
+    MESSAGE(STATUS "Could not find FFMPEG")
+
+ENDIF()
diff --git a/modules/FindLibMad.cmake b/modules/FindLibMad.cmake
new file mode 100644
index 0000000..20ad7d5
--- /dev/null
+++ b/modules/FindLibMad.cmake
@@ -0,0 +1,19 @@
+# Defines
+# LIBMAD_FOUND - if found
+# LIBMAD_INCLUDE_DIR - the include directory
+# LIBMAD_LIBRARY - libmad library path
+
+include(FindPackageHandleStandardArgs)
+
+find_path(LIBMAD_INCLUDE_DIR mad.h)
+find_library(LIBMAD_LIBRARY mad)
+
+find_package_handle_standard_args(
+  LibMad
+  DEFAULT_MSG
+  LIBMAD_LIBRARY
+  LIBMAD_INCLUDE_DIR
+)
+
+
+
diff --git a/paper/paper.aux b/paper/paper.aux
new file mode 100644
index 0000000..11a22f3
--- /dev/null
+++ b/paper/paper.aux
@@ -0,0 +1,6 @@
+\relax 
+\@writefile{toc}{\contentsline {section}{\numberline {1}Introduction}{1}}
+\@writefile{toc}{\contentsline {section}{\numberline {2}Tuning}{2}}
+\@writefile{toc}{\contentsline {section}{\numberline {3}Tempo Detection}{2}}
+\@writefile{toc}{\contentsline {section}{\numberline {4}Harmonic Analysis}{3}}
+\@writefile{toc}{\contentsline {section}{\numberline {5}Machine Learning}{3}}
diff --git a/paper/paper.log b/paper/paper.log
new file mode 100644
index 0000000..bbdfeb2
--- /dev/null
+++ b/paper/paper.log
@@ -0,0 +1,141 @@
+This is pdfTeX, Version 3.1415926-2.4-1.40.13 (TeX Live 2012) (format=pdflatex 2012.6.30)  4 AUG 2016 14:42
+entering extended mode
+ restricted \write18 enabled.
+ %&-line parsing enabled.
+**paper.tex
+(./paper.tex
+LaTeX2e <2011/06/27>
+Babel <v3.8m> and hyphenation patterns for english, dumylang, nohyphenation, ge
+rman-x-2012-05-30, ngerman-x-2012-05-30, afrikaans, ancientgreek, ibycus, arabi
+c, armenian, basque, bulgarian, catalan, pinyin, coptic, croatian, czech, danis
+h, dutch, ukenglish, usenglishmax, esperanto, estonian, ethiopic, farsi, finnis
+h, french, friulan, galician, german, ngerman, swissgerman, monogreek, greek, h
+ungarian, icelandic, assamese, bengali, gujarati, hindi, kannada, malayalam, ma
+rathi, oriya, panjabi, tamil, telugu, indonesian, interlingua, irish, italian, 
+kurmanji, latin, latvian, lithuanian, mongolian, mongolianlmc, bokmal, nynorsk,
+ polish, portuguese, romanian, romansh, russian, sanskrit, serbian, serbianc, s
+lovak, slovenian, spanish, swedish, turkish, turkmen, ukrainian, uppersorbian, 
+welsh, loaded.
+(/usr/local/texlive/2012/texmf-dist/tex/latex/base/article.cls
+Document Class: article 2007/10/19 v1.4h Standard LaTeX document class
+(/usr/local/texlive/2012/texmf-dist/tex/latex/base/size10.clo
+File: size10.clo 2007/10/19 v1.4h Standard LaTeX file (size option)
+)
+\c@part=\count79
+\c@section=\count80
+\c@subsection=\count81
+\c@subsubsection=\count82
+\c@paragraph=\count83
+\c@subparagraph=\count84
+\c@figure=\count85
+\c@table=\count86
+\abovecaptionskip=\skip41
+\belowcaptionskip=\skip42
+\bibindent=\dimen102
+) (./style.sty
+Package: style 
+)
+(/usr/local/texlive/2012/texmf-dist/tex/latex/amsmath/amsmath.sty
+Package: amsmath 2000/07/18 v2.13 AMS math features
+\@mathmargin=\skip43
+
+For additional information on amsmath, use the `?' option.
+(/usr/local/texlive/2012/texmf-dist/tex/latex/amsmath/amstext.sty
+Package: amstext 2000/06/29 v2.01
+
+(/usr/local/texlive/2012/texmf-dist/tex/latex/amsmath/amsgen.sty
+File: amsgen.sty 1999/11/30 v2.0
+\@emptytoks=\toks14
+\ex@=\dimen103
+))
+(/usr/local/texlive/2012/texmf-dist/tex/latex/amsmath/amsbsy.sty
+Package: amsbsy 1999/11/29 v1.2d
+\pmbraise@=\dimen104
+)
+(/usr/local/texlive/2012/texmf-dist/tex/latex/amsmath/amsopn.sty
+Package: amsopn 1999/12/14 v2.01 operator names
+)
+\inf@bad=\count87
+LaTeX Info: Redefining \frac on input line 211.
+\uproot@=\count88
+\leftroot@=\count89
+LaTeX Info: Redefining \overline on input line 307.
+\classnum@=\count90
+\DOTSCASE@=\count91
+LaTeX Info: Redefining \ldots on input line 379.
+LaTeX Info: Redefining \dots on input line 382.
+LaTeX Info: Redefining \cdots on input line 467.
+\Mathstrutbox@=\box26
+\strutbox@=\box27
+\big@size=\dimen105
+LaTeX Font Info:    Redeclaring font encoding OML on input line 567.
+LaTeX Font Info:    Redeclaring font encoding OMS on input line 568.
+\macc@depth=\count92
+\c@MaxMatrixCols=\count93
+\dotsspace@=\muskip10
+\c@parentequation=\count94
+\dspbrk@lvl=\count95
+\tag@help=\toks15
+\row@=\count96
+\column@=\count97
+\maxfields@=\count98
+\andhelp@=\toks16
+\eqnshift@=\dimen106
+\alignsep@=\dimen107
+\tagshift@=\dimen108
+\tagwidth@=\dimen109
+\totwidth@=\dimen110
+\lineht@=\dimen111
+\@envbody=\toks17
+\multlinegap=\skip44
+\multlinetaggap=\skip45
+\mathdisplay@stack=\toks18
+LaTeX Info: Redefining \[ on input line 2666.
+LaTeX Info: Redefining \] on input line 2667.
+) (./paper.aux)
+\openout1 = `paper.aux'.
+
+LaTeX Font Info:    Checking defaults for OML/cmm/m/it on input line 3.
+LaTeX Font Info:    ... okay on input line 3.
+LaTeX Font Info:    Checking defaults for T1/cmr/m/n on input line 3.
+LaTeX Font Info:    ... okay on input line 3.
+LaTeX Font Info:    Checking defaults for OT1/cmr/m/n on input line 3.
+LaTeX Font Info:    ... okay on input line 3.
+LaTeX Font Info:    Checking defaults for OMS/cmsy/m/n on input line 3.
+LaTeX Font Info:    ... okay on input line 3.
+LaTeX Font Info:    Checking defaults for OMX/cmex/m/n on input line 3.
+LaTeX Font Info:    ... okay on input line 3.
+LaTeX Font Info:    Checking defaults for U/cmr/m/n on input line 3.
+LaTeX Font Info:    ... okay on input line 3.
+ (./tex/title.tex) (./paper.toc)
+\tf@toc=\write3
+\openout3 = `paper.toc'.
+
+ (./tex/introduction.tex) (./tex/tuning.tex
+[1
+
+{/usr/local/texlive/2012/texmf-var/fonts/map/pdftex/updmap/pdftex.map}])
+(./tex/tempo.tex) (./tex/chords.tex [2]) (./tex/learning.tex) [3] (./paper.aux)
+ ) 
+Here is how much of TeX's memory you used:
+ 860 strings out of 493488
+ 9600 string characters out of 3141326
+ 63121 words of memory out of 3000000
+ 4199 multiletter control sequences out of 15000+200000
+ 8716 words of font info for 32 fonts, out of 3000000 for 9000
+ 957 hyphenation exceptions out of 8191
+ 32i,8n,22p,323b,189s stack positions out of 5000i,500n,10000p,200000b,50000s
+</usr/local/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmbx10.pfb
+></usr/local/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmbx12.pfb>
+</usr/local/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi10.pfb><
+/usr/local/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmr10.pfb></u
+sr/local/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmr12.pfb></usr
+/local/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmr17.pfb></usr/l
+ocal/texlive/2012/texmf-dist/fonts/type1/public/amsfonts/cm/cmr7.pfb>
+Output written on paper.pdf (3 pages, 90953 bytes).
+PDF statistics:
+ 42 PDF objects out of 1000 (max. 8388607)
+ 29 compressed objects within 1 object stream
+ 0 named destinations out of 1000 (max. 500000)
+ 1 words of extra memory for PDF output out of 10000 (max. 10000000)
+
diff --git a/paper/paper.pdf b/paper/paper.pdf
new file mode 100644
index 0000000..19e2c05
Binary files /dev/null and b/paper/paper.pdf differ
diff --git a/paper/paper.tex b/paper/paper.tex
new file mode 100644
index 0000000..ca5f9ae
--- /dev/null
+++ b/paper/paper.tex
@@ -0,0 +1,14 @@
+\documentclass{article}
+\usepackage{style}
+\begin{document}
+
+\input{./tex/title}
+\tableofcontents
+
+\input{./tex/introduction}
+\input{./tex/tuning}
+\input{./tex/tempo}
+\input{./tex/chords}
+\input{./tex/learning}
+
+\end{document}
diff --git a/paper/paper.toc b/paper/paper.toc
new file mode 100644
index 0000000..f3e40a6
--- /dev/null
+++ b/paper/paper.toc
@@ -0,0 +1,5 @@
+\contentsline {section}{\numberline {1}Introduction}{1}
+\contentsline {section}{\numberline {2}Tuning}{2}
+\contentsline {section}{\numberline {3}Tempo Detection}{2}
+\contentsline {section}{\numberline {4}Harmonic Analysis}{3}
+\contentsline {section}{\numberline {5}Machine Learning}{3}
diff --git a/paper/style.sty b/paper/style.sty
new file mode 100644
index 0000000..9619abb
--- /dev/null
+++ b/paper/style.sty
@@ -0,0 +1,3 @@
+\ProvidesPackage{style}
+
+\usepackage{amsmath}
diff --git a/paper/tex/chords.tex b/paper/tex/chords.tex
new file mode 100644
index 0000000..907d86e
--- /dev/null
+++ b/paper/tex/chords.tex
@@ -0,0 +1,7 @@
+\section{Harmonic Analysis}
+Using the tempo and the one as a guide,
+divide the song according to beat positions
+
+Equal loudness contour
+Normalizing energy per beat.
+Wrap to an octave (hard quantize, filter out very low signals)
diff --git a/paper/tex/introduction.tex b/paper/tex/introduction.tex
new file mode 100644
index 0000000..280934a
--- /dev/null
+++ b/paper/tex/introduction.tex
@@ -0,0 +1,18 @@
+\section{Introduction}
+
+Sample based music has an 
+Samples of Motown backing bands, pop stars, and jazz legends can be stitched together to create a collage of sound combining all of their talents.
+Sampling has the potential to make music that was never possible before.
+Out of a musician's life work, the sampler can choose a 10 second snippet that perfectly
+Condensing hours and hours of recording time, studios across the globe and talents of all sorts.
+
+While some artists have made sampled based masterpieces, they take considerable ammount of time.
+The issue
+The goal of this project is to reduce almost all of the manual labor associated with sampling.
+The process should ideally be simply listening to music, choosing what sounds good, and then
+In reality lots of time is spent tuning samples, matching beats, and matching chords.
+
+Phase vocording works to solve this problem, it nessesarily introduces digital artifacts.
+Although subtle, these artifacts can break the suspension of disbelief holding together the dream that Madonna and James Brown really did record together.
+Rather than using computers to directly manipulate audio, this work intends to do an intense analysis on a sample library, beyond the scope of humans, to find pairings for samples.
+While more digital processing has been done to the audio than probably any other production, the actual audio signal path could be completely analog. 
diff --git a/paper/tex/learning.tex b/paper/tex/learning.tex
new file mode 100644
index 0000000..5fc2eac
--- /dev/null
+++ b/paper/tex/learning.tex
@@ -0,0 +1,14 @@
+\section{Machine Learning}
+
+ongs that work together
+Scale polynomials (tune) and find starting point,
+so that fluctation of beat over length is negligible
+
+Recurrent neural nets and PU Learning
+Train the NN to learn the difference between positive examples
+(real songs)
+and songs that have been built like above
+The samples should be
+Punish false positives less than false negatives.
+So that all positive examples are learned, and some negative ones will bleed through
+Search for ones that do indeed bleed through.
diff --git a/paper/tex/tempo.tex b/paper/tex/tempo.tex
new file mode 100644
index 0000000..19c6be7
--- /dev/null
+++ b/paper/tex/tempo.tex
@@ -0,0 +1,69 @@
+\section{Tempo Detection}
+
+Tempo detection
+Approximating the tempo
+Onset detection
+Autocorrelation
+Then use wavelets to detect self similarity
+
+Tempo detection happens in three steps.
+First we convert the audio signal into an ``onset'' signal, which is a positive signal whose amplitude correlates to to how much the signal is changing at that point. This signal is high when the amplitude or pitch of the signal suddenly changes, which would be a ``beat.''
+Next we approximate the tempo by looking for the frequency of self simiarity in this onset signal.
+Finally, we fine tune this approximation to get function for the tempo that minimizes the ammount of onset information that occurs off-beat.
+
+Onset signal
+
+Using an analysis by [], we are going to 
+We window the audio signal into frames.
+We use the previous two frames to linearly predict what the third frame should be.
+If the signal is changing this third frame will be different than the previous two.
+So we simply calculate the euclidean distance between our predicted frame and the actual frame.
+
+Approximating the signal
+
+A tempo in music is the frequency at which the music exhibits some sort of self similarity. Self similariy on music exists on many levels. There is repeating verse structure - 16 bars, repeating chord progressions - 4 bars, repeating drum beats - 1 bar, and the beat itself - 1 quarter note.
+In order to divide up a song as finely as to display its harmonic and rythmic features want to find the minimum frequency at which the song exhibits self similarity. 
+
+Fine tuning the approximation
+
+Many modern music recordings are 
+
+While in most cases playing to a click or quantized electronic music has a perfectly consistent tempo, we would like to include music played live, where the tempo can slowly shift or change
+Series approximation with wavelets?
+Fourier series?
+Polynomial approximation
+Approximation
+$b(t) = a_0+a_1t+...$
+$b(t) = a_0 + a_1cos($
+Find coeffiecient, then next coefficient must have average of first ...
+Most songs will have 0 higher order coefficients
+But those that fluctuate will not.
+
+We want to minimize the function
+\begin{align*}
+  F =
+  \sum_{bin = 0}^{numBins}
+  \eta(bin)(1 - \cos(f(t))
+\end{align*}
+Where
+\begin{align*}
+  f(t) = \alpha t + a_0 + 
+  \sum_{i = 1}^k 
+  a_i\cos\left(\frac{2i\pi t}{N}\right)
+  +
+  b_i\sin\left(\frac{2i\pi t}{N}\right)
+\end{align*}
+
+\begin{align*}
+  \nabla F = 
+  \sum_{bin = 0}^{numBins}
+  \left(
+  -\eta(bin)\sin(f(t))
+  +
+  \left( 
+  \sum_{i = 1}^k
+  \eta(bin)\frac{N}{2i\pi t}\sin\left(\frac{2i\pi t}{N}\right)\sin(f(t))
+  -\eta(bin)\frac{N}{2i\pi t}\cos\left(\frac{2i\pi t}{N}\right)\sin(f(t))
+  \right)
+  \right)
+\end{align*}
diff --git a/paper/tex/title.tex b/paper/tex/title.tex
new file mode 100644
index 0000000..47ee611
--- /dev/null
+++ b/paper/tex/title.tex
@@ -0,0 +1,4 @@
+\title{Sampling}
+\date{\today}
+\author{Trevor Henderson}
+\maketitle
diff --git a/paper/tex/tuning.tex b/paper/tex/tuning.tex
new file mode 100644
index 0000000..e682ac5
--- /dev/null
+++ b/paper/tex/tuning.tex
@@ -0,0 +1,14 @@
+\section{Tuning}
+
+In order to get the sample ready to be heard by our machine learner, the first thing we do is tune the sample so that the frequencies will fit nicely into semitone bins.
+This is the easiest step, but it introduces some techniques we will be using later on.
+
+We will assume, very reasonably, that the song being taken as the input is not long enough that the instruments began to go out of tune, and was recorded and transcoded with good technology - a turntable with a stable belt, etc. 
+Therefore the tuning will be constant which makes our job easier.
+
+What we will do is take the fourier transform and then lump frequencies together into bins. 
+This is fairly good for large samples, but for small samples, sometimes the frequency resolution isn't fine enough to discriminate between semitones, so we will use frequency estimation:
+
+\begin{align*}
+  \kappa(x) = 
+\end{align*}
diff --git a/paper/tex/writing.tex b/paper/tex/writing.tex
new file mode 100644
index 0000000..8601a8e
--- /dev/null
+++ b/paper/tex/writing.tex
@@ -0,0 +1,86 @@
+Musical analysis
+
+Sample based music has an 
+Samples of Motown backing bands, pop stars, and jazz legends can be stitched together to create a collage of sound combining all of their talents.
+Sampling has the potential to make music that was never possible before.
+Out of a musician's life work, the sampler can choose a 10 second snippet that perfectly
+Condensing hours and hours of recording time, studios across the globe and talents of all sorts.
+
+While some artists have made sampled based masterpieces, they take considerable ammount of time.
+The issue
+The goal of this project is to reduce almost all of the manual labor associated with sampling.
+The process should ideally be simply listening to music, choosing what sounds good, and then
+In reality lots of time is spent tuning samples, matching beats, and matching chords.
+
+Phase vocording works to solve this problem, it nessesarily introduces digital artifacts.
+Although subtle, these artifacts can break the suspension of disbelief holding together the dream that Madonna and James Brown really did record together.
+Rather than using computers to directly manipulate audio, this work intends to do an intense analysis on a sample library, beyond the scope of humans, to find pairings for samples.
+While more digital processing has been done to the audio than probably any other production, the actual audio signal path could be completely analog. 
+
+
+Tuning
+This is a constant. We assume the song is short enough and was transcoded though suitable means so it does not detune while it is being played. 
+
+Frequency estimation
+amplitude estimation
+
+Tempo detection
+Approximating the tempo
+Onset detection
+Autocorrelation
+Then use wavelets to detect self similarity
+
+Tempo detection happens in three steps.
+First we convert the audio signal into an ``onset'' signal, which is a positive signal whose amplitude correlates to to how much the signal is changing at that point. This signal is high when the amplitude or pitch of the signal suddenly changes, which would be a ``beat.''
+Next we approximate the tempo by looking for the frequency of self simiarity in this onset signal.
+Finally, we fine tune this approximation to get function for the tempo that minimizes the ammount of onset information that occurs off-beat.
+
+Onset signal
+
+Using an analysis by [], we are going to 
+We window the audio signal into frames.
+We use the previous two frames to linearly predict what the third frame should be.
+If the signal is changing this third frame will be different than the previous two.
+So we simply calculate the euclidean distance between our predicted frame and the actual frame.
+
+Approximating the signal
+
+A tempo in music is the frequency at which the music exhibits some sort of self similarity. Self similariy on music exists on many levels. There is repeating verse structure - 16 bars, repeating chord progressions - 4 bars, repeating drum beats - 1 bar, and the beat itself - 1 quarter note.
+In order to divide up a song as finely as to display its harmonic and rythmic features want to find the minimum frequency at which the song exhibits self similarity. 
+
+Fine tuning the approximation
+
+Many modern music recordings are 
+
+While in most cases playing to a click or quantized electronic music has a perfectly consistent tempo, we would like to include music played live, where the tempo can slowly shift or change
+Series approximation with wavelets?
+Fourier series?
+Polynomial approximation
+Approximation
+b(t) = a_0+a_1t+...
+b(t) = a_0 + a_1cos(
+Find coeffiecient, then next coefficient must have average of first ...
+Most songs will have 0 higher order coefficients
+But those that fluctuate will not.
+
+Chords
+Using the tempo and the one as a guide,
+divide the song according to beat positions
+
+Equal loudness contour
+Normalizing energy per beat.
+Wrap to an octave (hard quantize, filter out very low signals)
+
+Finding songs that work together
+Scale polynomials (tune) and find starting point,
+so that fluctation of beat over length is negligible
+
+Recurrent neural nets and PU Learning
+Train the NN to learn the difference between positive examples
+(real songs)
+and songs that have been built like above
+The samples should be
+Punish false positives less than false negatives.
+So that all positive examples are learned, and some negative ones will bleed through
+Search for ones that do indeed bleed through.
+
diff --git a/python_src/ProcessLibrary.py b/python_src/ProcessLibrary.py
index b5a78b8..f1615e7 100644
--- a/python_src/ProcessLibrary.py
+++ b/python_src/ProcessLibrary.py
@@ -1,45 +1,157 @@
-import os
+import os, math, signal
 
 from sample import Sample
 
-library = "/Users/tfh/Dropbox (MIT)/UserLibrary/SampleLibrary/Chops/"
+#import tensorflow as tf
 
-samples = []
+#def getRandomChords(samples):
 
-for root, directories, filenames in os.walk(library):
-    for filename in filenames:
-        if filename[-4:] == ".alc":
-            filePath = os.path.join(root, filename)
-            print "Processing:", filename[:-4]
-            s = Sample(filePath)
-            sample.append(s)
-            print "Tuning:", s.getTuning()
-            print "Tempo:", s.getTempo() * 60.
-            print
+#def getRandomCHords(pairs):
 
-ratios = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 3/2.]
+# to generate a batch:
+# for i in batch size:
+# choose randomly between regular (1) and combined (-1)
+# if randomly choose the intersection point
+# randomly rotate the chords
+# add and average energy...
+# if combined randomly choose start points (its recurrent so train until one of them ends)
+# train the output to be 1 consistent
+# choose start point randomly and train to the end.
 
-for s1 in samples:
-    tempo1 = s1.getTempo()
-    for s2 in samples:
-        tempo2 = s2.getTempo()
+# listen to midi to see if euclidean distance is better than
+# simply adding up...
 
-        larger = max(tempo1, tempo2)
-        smaller = min(tempo1, tempo2)
+#lstm = tf.nn.rnn_cell.BasicLSTMCELL(lstm_size)
+#fw_cell = tf.nn.rnn_cell.MultiRNNCell([lstm] * number_of_layers)
+#bw_cell = fw_cell
 
-        ratio = large/small
+#initial_state = state = stacked_lstm.zero_state(None, tf.float32)
 
-        for r in ratios:
-            # if pitch change < an octave its cool
-             
 
 
+#x = Placeholder for list of tensors of size [None, 12]
 
-        # round to 2, 3, 4, 6, 8, 12, 16, 24, 32 or 3/2, 4/3, 9/2, 9/4
+#tf.nn.bidirectional_rnn(fw_cell, bw_cell, x)
 
-        # get larger of the two
-        # what is whole number ratio large/small
-        # how far away is it from that ratio
-        # also take into account ratios like 2/3
-        
-        
+def processFiles(library):
+    processedFiles = []
+    unProcessedFiles = []
+    for root, directories, filenames in os.walk(library):
+        for filename in filenames:
+              if filename[-4:] == ".alc":
+                  filePath = os.path.join(root, filename)
+                  s = Sample(filePath)
+                  if s.process():
+                      print "Processed:", s.getFileName()
+                      print "Tuning:", s.getTuning()
+                      print "Tempo:", s.getTempo() * 60.
+                      print
+                      s.writeToFile()
+                      s.delete()
+                      processedFiles.append(filePath)
+                  else:
+                      print "Could not process ", filename
+                      print
+    return processedFiles
+
+def findPairs(files, semitonesLowerBound, semitonesUpperBound, tuningBound):
+    ratios = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]
+    inverseRatios = [1./r for r in ratios]
+    ratios += inverseRatios
+
+    print ratios
+
+    pairs = {}
+
+    for i in range(len(files)):
+        s1 = Sample(files[i])
+        s1.process()
+        tempo1 = s1.getTempo()
+
+        pairs[s1.getFileName()] = []
+
+        for j in range(i+1, len(files)):
+            s2 = Sample(files[j])
+            s2.process()
+            tempo2 = s2.getTempo()
+
+            ratio = tempo1/tempo2
+
+            for r in ratios:
+                numCents = 1200. * math.log(ratio/r, 2.)
+                distanceFromSemitone = (numCents % 100.)
+                if distanceFromSemitone > 50:
+                    distanceFromSemitone = 100 - distanceFromSemitone
+
+                if semitonesLowerBound * 100. < numCents < semitonesUpperBound*100. \
+                    and distanceFromSemitone < tuningBound:
+                    pairs[s1.getFileName()].append((s2.getFileName(), r, numCents))
+            s2.delete()
+        s1.delete()
+
+    return pairs
+
+def makeBatch(batchSize, files, pairs):
+  batch = []
+  for i in range(batchSize):
+    # choose randomly between file of pair
+    # if file:
+    # choose a random file 
+    # get chords as vector
+    # if pair
+    # choose a random pair
+    # get chords for both
+    # choose a random intersection point
+    # add them up and divide by 2
+    # for both:
+    # randomly rotate
+    pass
+  return batch
+
+
+def main():
+    library = "/Users/tfh/Dropbox (MIT)/UserLibrary/SampleLibrary/Chops/"
+    #library = "/Users/tfh/Dropbox (MIT)/SampleSorter/TestFiles/"
+
+    processedFiles = processFiles(library)
+
+    pairs = findPairs(processedFiles, -6, 3, 5)
+
+    print len(pairs),"pairs!"
+
+    maxPairs = 0
+    for f in pairs:
+      pairsF = len(pairs[f])
+      if pairsF > maxPairs:
+        maxPairs = pairsF
+        maxF = f
+
+    for fileName, ratio, numCents in pairs[maxF]:
+      sample = Sample(fileName)
+      print sample.getFileName()
+      print ratio
+      print numCents
+      sample.delete()
+
+    samplF = Sample(maxF)
+    print samplF.getFileName(), "has", maxPairs,"pairs!"
+    samplF.delete()
+
+    
+
+
+# TODO
+# GUI
+# find cliques
+# drag and drop
+# drop files in
+# drag files out with appropriate relative tuning
+# mp3 support
+# warped files support
+# straight audio file support
+# support for the couple samples that aren't working
+
+
+
+if __name__ == "__main__":
+  main()
diff --git a/python_src/avalanche2.mid b/python_src/avalanche2.mid
new file mode 100644
index 0000000..9055a08
Binary files /dev/null and b/python_src/avalanche2.mid differ
diff --git a/python_src/bach.mid b/python_src/bach.mid
index e4694e8..afca5af 100644
Binary files a/python_src/bach.mid and b/python_src/bach.mid differ
diff --git a/python_src/beachBoys.mid b/python_src/beachBoys.mid
new file mode 100644
index 0000000..5f07802
Binary files /dev/null and b/python_src/beachBoys.mid differ
diff --git a/python_src/joe3.mid b/python_src/joe3.mid
new file mode 100644
index 0000000..be397f1
Binary files /dev/null and b/python_src/joe3.mid differ
diff --git a/python_src/sample.py b/python_src/sample.py
index 10b416a..00e68ba 100644
--- a/python_src/sample.py
+++ b/python_src/sample.py
@@ -9,6 +9,16 @@ class Sample:
         lib.NewAbletonSampleFile.restype = ctypes.c_void_p
         self.s = lib.NewAbletonSampleFile(str(fileName).encode('ascii'))
 
+    def getFileName(self):
+        lib.getFileName.argtypes = [ctypes.c_void_p]
+        lib.getFileName.restype = ctypes.c_char_p
+        return lib.getFileName(self.s)
+        
+    def process(self):
+        lib.process.argtypes = [ctypes.c_void_p]
+        lib.process.restype = ctypes.c_bool
+        return lib.process(self.s)
+
     def getTuning(self):
         lib.getTuningCents.argtypes = [ctypes.c_void_p]
         lib.getTuningCents.restype = ctypes.c_long
@@ -43,6 +53,10 @@ class Sample:
         # return
         return chords
 
+    def writeToFile(self):
+        lib.writeToFile.argtypes = [ctypes.c_void_p]
+        lib.writeToFile(self.s)
+
     def delete(self):
         lib.deleteAbletonSampleFile.argtypes = [ctypes.c_void_p]
         lib.deleteAbletonSampleFile(self.s)
diff --git a/python_src/sample.pyc b/python_src/sample.pyc
index e4cce86..2ef08fa 100644
Binary files a/python_src/sample.pyc and b/python_src/sample.pyc differ
diff --git a/src/AbletonSampleFile.cpp b/src/AbletonSampleFile.cpp
index dd0ffe9..a8b6e50 100644
--- a/src/AbletonSampleFile.cpp
+++ b/src/AbletonSampleFile.cpp
@@ -1,8 +1,10 @@
 #include <iostream>
 
 #include <vector>
+#include <string>
 #include <fstream>
 #include <sstream>
+#include <initializer_list>
 
 #include <boost/iostreams/filter/gzip.hpp>
 #include <boost/iostreams/filtering_stream.hpp>
@@ -11,10 +13,10 @@
 
 #include <tinyxml2.h>
 
-#include <sndfile.hh>
-
 #include "SampleSorter/SampleFile.hpp"
 #include "SampleSorter/AbletonSampleFile.hpp"
+#include "SampleSorter/ProcessingException.hpp"
+#include "SampleSorter/AudioFile.hpp"
 
 AbletonSampleFile::AbletonSampleFile(const AbletonSampleFile & other) 
   : SampleFile(other)
@@ -24,8 +26,11 @@ AbletonSampleFile::AbletonSampleFile(const AbletonSampleFile & other)
 
 AbletonSampleFile::AbletonSampleFile(std::string filePath) 
   : SampleFile(filePath) {
-    getDoc();
-    readDoc();
+}
+
+bool AbletonSampleFile::readMetaData() {
+  getDoc();
+  return readDoc();
 }
 
 std::string AbletonSampleFile::getReferenceFilePath() const {
@@ -36,97 +41,273 @@ std::string AbletonSampleFile::getReferenceFileName() const {
   return boost::filesystem::path(referenceFilePath).stem().native();
 }
 
-double AbletonSampleFile::getSampleLength() const {
+double AbletonSampleFile::getSampleSeconds() const {
   return endSeconds - startSeconds;
 }
 
-long AbletonSampleFile::getSampleRate() const {
-  return audioFile.samplerate();
+std::vector< std::vector<double> > AbletonSampleFile::extractAudio(long * sampleRate) {
+  return AudioFile::read(referenceFilePath, startSeconds, endSeconds, sampleRate);
+}
+
+void AbletonSampleFile::getDoc() {
+  std::stringstream unzipped;
+  try {
+    // Ungzip
+    std::ifstream zipped(filePath, 
+                         std::ios_base::in | std::ios_base::binary);
+    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
+    in.push(boost::iostreams::gzip_decompressor());
+    in.push(zipped);
+    boost::iostreams::copy(in, unzipped);
+  } catch (boost::iostreams::gzip_error & e) {
+    throw ProcessingException("Could not ungzip Ableton file!");
+  }
+
+  // convert from XML
+  doc.Parse(&unzipped.str()[0]);
 }
 
-std::vector< std::vector<double> > AbletonSampleFile::getWaves() {
-  // set seek to beginning
-  audioFile.seek(startSeconds * getSampleRate(), 0);
-  // number of frames to read
-  long size = (endSeconds - startSeconds) * getSampleRate();
-  std::vector<double> rawAudioData(size * audioFile.channels());
-  // read raw interleaved channels
-  audioFile.read(&rawAudioData[0], size * audioFile.channels());
+void AbletonSampleFile::setDoc() {
+  // make a string
+  tinyxml2::XMLPrinter printer;
+  doc.Accept(&printer);
 
-  // Allocate output vector
-  std::vector<std::vector<double> > waves(audioFile.channels());
-  for (int channel = 0; channel < audioFile.channels(); channel ++) {
-    waves[channel].resize(size);
+  // gzip
+  try {
+    std::stringstream unzipped;
+    unzipped << printer.CStr();
+
+    std::ofstream zipped(filePath, 
+                         std::ios_base::out | std::ios_base::binary);
+    boost::iostreams::filtering_streambuf<boost::iostreams::input> filter;
+    filter.push(boost::iostreams::gzip_compressor());
+    filter.push(unzipped);
+    boost::iostreams::copy(filter, zipped);
+  } catch (boost::iostreams::gzip_error & e) {
+    throw ProcessingException("Could not gzip Ableton file!");
   }
+}
 
-  // De-interleave channels
-  for (long i = 0; i < size * audioFile.channels(); i++) {
-    waves[i % audioFile.channels()][i/audioFile.channels()] = rawAudioData[i];
+tinyxml2::XMLElement * AbletonSampleFile::getNode(
+    std::initializer_list<std::string> nodes
+    ) {
+  tinyxml2::XMLElement * node = doc.RootElement();
+  for (auto it = nodes.begin(); it < nodes.end(); it++) {
+    node = node -> FirstChildElement((*it).c_str());
   }
+  return node;
+}
 
-  return waves;
+tinyxml2::XMLElement * AbletonSampleFile::getAudioNode() {
+  return getNode({
+      "LiveSet",
+      "Tracks",
+      "AudioTrack",
+      "DeviceChain",
+      "MainSequencer",
+      "ClipSlotList",
+      "ClipSlot",
+      "ClipSlot",
+      "Value",
+      "AudioClip"
+      });
 }
 
-void AbletonSampleFile::getDoc() {
-  // Ungzip
-  std::stringstream unzipped;
+tinyxml2::XMLElement * AbletonSampleFile::getLoopNode() {
+  return getAudioNode() -> FirstChildElement("Loop");
+}
 
-  std::ifstream zipped(filePath, 
-                       std::ios_base::in | std::ios_base::binary);
-  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
-  in.push(boost::iostreams::gzip_decompressor());
-  in.push(zipped);
-  boost::iostreams::copy(in, unzipped);
+tinyxml2::XMLElement * AbletonSampleFile::getSortDataNode() {
+  return getNode({"LiveSet","SampleSorter"});
+}
 
-  // convert from XML
-  doc.Parse(&unzipped.str()[0]);
+tinyxml2::XMLElement * AbletonSampleFile::getLoopStartNode() {
+  return getLoopNode() -> FirstChildElement("LoopStart");
+}
+
+tinyxml2::XMLElement * AbletonSampleFile::getLoopEndNode() {
+  return getLoopNode() -> FirstChildElement("LoopEnd");
 }
 
-void AbletonSampleFile::readDoc() {
-  tinyxml2::XMLElement * audioNode = doc.RootElement()
-                                          -> FirstChildElement("LiveSet")
-                                          -> FirstChildElement("Tracks")
-                                          -> FirstChildElement("AudioTrack")
-                                          -> FirstChildElement("DeviceChain")
-                                          -> FirstChildElement("MainSequencer")
-                                          -> FirstChildElement("ClipSlotList")
-                                          -> FirstChildElement("ClipSlot")
-                                          -> FirstChildElement("ClipSlot")
-                                          -> FirstChildElement("Value")
-                                          -> FirstChildElement("AudioClip");
-  tinyxml2::XMLElement * loopNode = audioNode -> FirstChildElement("Loop");
+tinyxml2::XMLElement * AbletonSampleFile::getHiddenLoopStartNode() {
+  return getLoopNode() -> FirstChildElement("HiddenLoopStart");
+}
 
-  // Start and end times in seconds
-  startSeconds = loopNode -> FirstChildElement("LoopStart") 
-                                -> DoubleAttribute("Value");
-  endSeconds = loopNode -> FirstChildElement("LoopEnd") 
-                                -> DoubleAttribute("Value");
-
-  tinyxml2::XMLElement * pathNode = audioNode -> FirstChildElement("SampleRef")
-                                              -> FirstChildElement("FileRef");
-  if (pathNode == nullptr) {
-    pathNode = audioNode -> FirstChildElement("SampleRef")
-                         -> FirstChildElement("SourceContext")
-                         -> FirstChildElement("SourceContext")
-                         -> FirstChildElement("OriginalFileRef")
-                         -> FirstChildElement("FileRef");
+tinyxml2::XMLElement * AbletonSampleFile::getHiddenLoopEndNode() {
+  return getLoopNode() -> FirstChildElement("HiddenLoopEnd");
+}
+
+tinyxml2::XMLElement * AbletonSampleFile::getNameNode() {
+  return getAudioNode() -> FirstChildElement("Name");
+}
+
+tinyxml2::XMLElement * AbletonSampleFile::getPitchFineNode() {
+  return getAudioNode() -> FirstChildElement("PitchFine");
+}
+
+std::string AbletonSampleFile::fetchReferenceFilePath() {
+  // Some awful reverse engineering
+  std::string path;
+
+  tinyxml2::XMLElement * pathNode = getAudioNode() 
+    -> FirstChildElement("SampleRef")
+    -> FirstChildElement("SourceContext")
+    -> FirstChildElement("SourceContext");
+
+  if (pathNode != nullptr) {
+    pathNode = pathNode -> FirstChildElement("BrowserContentPath");
+
+    path = pathNode -> Attribute("Value");
+    // remove "userlibrary:"
+    path.erase(0,11);
+
+    // replace "%20" with spaces
+    size_t replaceIndex = 0;
+    while (true) {
+      replaceIndex = path.find("%20", replaceIndex);
+      if (replaceIndex == std::string::npos) break;
+
+      path.replace(replaceIndex, 3, " ");
+      replaceIndex += 1;
+    }
+    replaceIndex = 0;
+
+    // remove the first #
+    replaceIndex = path.find("#", replaceIndex);
+    path.erase(replaceIndex, 1);
+
+    // replace the :s with /
+    replaceIndex = 0;
+    while (true) {
+      replaceIndex = path.find(":", replaceIndex);
+      if (replaceIndex == std::string::npos) break;
+
+      path.replace(replaceIndex, 1, "/");
+      replaceIndex += 1;
+    }
+  } else { 
+    pathNode = getAudioNode() -> FirstChildElement("SampleRef")
+                              -> FirstChildElement("FileRef");
+
+    tinyxml2::XMLElement * dirNode = pathNode -> FirstChildElement("SearchHint")
+                                      -> FirstChildElement("PathHint")
+                                      -> FirstChildElement("RelativePathElement");
+
+    // get file path
+    path = "/";
+    while (dirNode != nullptr) {
+      path += dirNode -> Attribute("Dir");
+      path += "/";
+      dirNode = dirNode -> NextSiblingElement("RelativePathElement");
+    }
+
+    // get file name
+    tinyxml2::XMLElement * nameNode = pathNode -> FirstChildElement("Name");
+    path += nameNode -> Attribute("Value");
   }
 
-  tinyxml2::XMLElement * dirNode = pathNode -> FirstChildElement("SearchHint")
-                                    -> FirstChildElement("PathHint")
-                                    -> FirstChildElement("RelativePathElement");
+  return path;
+}
 
-  // get file path
-  referenceFilePath = "/";
-  while (dirNode != nullptr) {
-    referenceFilePath += dirNode -> Attribute("Dir");
-    referenceFilePath += "/";
-    dirNode = dirNode -> NextSiblingElement("RelativePathElement");
+bool AbletonSampleFile::readDoc() {
+
+  // Start and end times in seconds
+  startSeconds = getLoopStartNode() -> DoubleAttribute("Value");
+  endSeconds = getLoopEndNode() -> DoubleAttribute("Value");
+
+  // if this has already been processed
+  if (getSortDataNode() != nullptr) {
+    // Path
+    referenceFilePath = getSortDataNode() -> Attribute("ReferenceFilePath");
+    // Tuning
+    long tuningCents = getPitchFineNode() -> IntAttribute("Value");
+    // Tempo
+    double rawBeat = getSortDataNode() -> DoubleAttribute("RawBeat");
+    // One
+    double theOne = getSortDataNode() -> DoubleAttribute("TheOne");
+    // SampleRate
+    double sampleRate = getSortDataNode() -> IntAttribute("SampleRate");
+
+    // Chords
+    std::vector<Octave> chords;
+    tinyxml2::XMLElement * chordElement = getSortDataNode() 
+      -> FirstChildElement("Chord");
+    while (chordElement != nullptr) {
+      std::vector<double> spectrogram(12);
+      tinyxml2::XMLElement * noteElement = 
+        chordElement -> FirstChildElement("Note");
+      for (long i = 0; i < 12; i++) {
+        spectrogram[i] = noteElement -> DoubleAttribute("Value");
+        noteElement = noteElement -> NextSiblingElement("Note");
+      }
+      chordElement = chordElement -> NextSiblingElement("Chord");
+      Octave newChord(spectrogram);
+      chords.push_back(newChord);
+    }
+
+    // Make the sample!
+    sample = AudioSample(
+        tuningCents, 
+        rawBeat, 
+        theOne, 
+        endSeconds - startSeconds, 
+        sampleRate,
+        chords
+        );
+
+    return true;
+  } else {
+    // it was not already processed, return false
+    referenceFilePath = fetchReferenceFilePath();
+    return false;
   }
+}
 
-  // get file name
-  tinyxml2::XMLElement * nameNode = pathNode -> FirstChildElement("Name");
-  referenceFilePath += nameNode -> Attribute("Value");
+void AbletonSampleFile::writeToFile() {
+  // Delete old data to write new
+  tinyxml2::XMLElement * sortData = getSortDataNode();
+  doc.RootElement() -> FirstChildElement("LiveSet")
+                    -> DeleteChild(sortData);
+  sortData = doc.NewElement("SampleSorter");
+
+  // Path
+  sortData -> SetAttribute("ReferenceFilePath", referenceFilePath.c_str());
+  // Tempo
+  sortData -> SetAttribute("RawBeat", getAudioSample() -> getBeatRaw());
+  // The one
+  sortData -> SetAttribute("TheOne", getAudioSample() -> getTheOneRaw());
+  // Sample rate
+  sortData -> SetAttribute("SampleRate", int(getAudioSample() -> getSampleRate()));
+
+  // Loop ends
+  getHiddenLoopStartNode() -> SetAttribute("Value", startSeconds);
+  getHiddenLoopEndNode() -> SetAttribute("Value", endSeconds);
+
+  // Tuning
+  getPitchFineNode() 
+    -> SetAttribute("Value", int(getAudioSample() -> getTuningCents()));
+
+  // Name
+  getNameNode() -> SetAttribute("Value",
+      (
+        getFileName() + " @" +
+        std::to_string(60. * (getAudioSample() -> getBeatWithTuning()))
+      ).c_str()
+    );
+
+  // Chords
+  std::vector<Octave> chords = getAudioSample() -> getChords();
+  for (long i = 0; i < chords.size(); i++) {
+    tinyxml2::XMLElement * chordElement = doc.NewElement("Chord");
+    for (short j = 0; j < 12; j++) {
+      tinyxml2::XMLElement * noteElement = doc.NewElement("Note");
+      noteElement -> SetAttribute("Value", chords[i].getSpectrogram()[j]);
+      chordElement -> InsertEndChild(noteElement);
+    }
+    sortData -> InsertEndChild(chordElement);
+  }
 
-  audioFile = SndfileHandle(referenceFilePath);
+  // Put in doc
+  doc.RootElement() -> FirstChildElement("LiveSet") -> InsertEndChild(sortData); 
+  setDoc();
 }
diff --git a/src/AudioFile.cpp b/src/AudioFile.cpp
new file mode 100644
index 0000000..96c8231
--- /dev/null
+++ b/src/AudioFile.cpp
@@ -0,0 +1,187 @@
+#include <iostream>
+
+#include <vector>
+#include <string>
+
+extern "C" {
+#include <libavformat/avformat.h>
+};
+
+#include "SampleSorter/AudioFile.hpp"
+#include "SampleSorter/ProcessingException.hpp"
+
+std::vector<std::vector<double> > AudioFile::read(
+    std::string fileName,
+    double startSeconds,
+    double endSeconds,
+    long * sampleRate) {
+
+  if (startSeconds < 0 or endSeconds < 0) {
+    throw ProcessingException("Time is negative!");
+  }
+
+  if (startSeconds > endSeconds) {
+    throw ProcessingException("The time to start reading is after the time to end reading");
+  }
+
+  // Initialize FFmpeg
+  av_register_all();
+
+  // Allocate a frame
+  AVFrame * frame = av_frame_alloc();
+  if (!frame) {
+    throw ProcessingException("Could not allocate audio frame!");
+  }
+
+  // open file
+  AVFormatContext * formatContext = NULL;
+  if (avformat_open_input(&formatContext, fileName.c_str(), NULL, NULL) != 0) {
+    av_free(frame);
+    throw ProcessingException("Could not open audio file!");
+  }
+
+  // get stream info
+  if (avformat_find_stream_info(formatContext, NULL) < 0) {
+    av_free(frame);
+    avformat_close_input(&formatContext);
+    throw ProcessingException("Could not get stream info!");
+  }
+
+  // Get the audio stream
+  AVCodec * cdc = nullptr;
+  int streamIndex = av_find_best_stream(formatContext, AVMEDIA_TYPE_AUDIO, -1, -1, &cdc, 0);
+  if (streamIndex < 0) {
+    // Could not find audio stream in file
+    av_free(frame);
+    avformat_close_input(&formatContext);
+    throw ProcessingException("Could not get stream info!");
+  }
+
+  AVStream * audioStream = formatContext -> streams[streamIndex];
+  AVCodecContext * codecContext = audioStream -> codec;
+  codecContext -> codec = cdc;
+
+  if (avcodec_open2(codecContext, codecContext -> codec, NULL) != 0) {
+    av_free(frame);
+    avformat_close_input(&formatContext);
+    throw ProcessingException("Could not open the context with decoder!");
+  }
+
+  *sampleRate = codecContext -> sample_rate;
+
+  long startSample = startSeconds * (*sampleRate);
+  long endSample = endSeconds * (*sampleRate);
+
+  std::vector<std::vector<double> > output(codecContext -> channels);
+  for (long channel = 0; channel < codecContext -> channels; channel++) {
+    output[channel].resize(endSample - startSample);
+  }
+
+  long sampleOffset = 0;
+
+  AVPacket readingPacket;
+  av_init_packet(&readingPacket);
+
+  while (sampleOffset < endSample and
+      av_read_frame(formatContext, &readingPacket) == 0) {
+
+    if (readingPacket.stream_index ==audioStream -> index) {
+      AVPacket decodingPacket = readingPacket;
+
+      while (decodingPacket.size > 0) {
+        int gotFrame = 0;
+        int result = avcodec_decode_audio4(codecContext, frame, &gotFrame, &decodingPacket);
+
+        if (result >= 0 && gotFrame) {
+          decodingPacket.size -= result;
+          decodingPacket.data += result;
+
+          if (sampleOffset + frame -> nb_samples > startSample) {
+            readFrame(codecContext, 
+                      frame, 
+                      output, 
+                      startSample, 
+                      endSample, 
+                      sampleOffset);
+          }
+
+          sampleOffset += frame -> nb_samples;
+
+        } else {
+          decodingPacket.size = 0;
+          decodingPacket.data = nullptr;
+        }
+      }
+    }
+    av_packet_unref(&readingPacket);
+  }
+
+  // Clean up
+  av_free(frame);
+  avcodec_close(codecContext);
+  avformat_close_input(&formatContext);
+
+  return output;
+}
+
+
+void AudioFile::readFrame(
+    const AVCodecContext * codecContext, 
+    const AVFrame * frame,
+    std::vector<std::vector<double> > & output,
+    long startSample,
+    long endSample,
+    long sampleOffset) {
+
+  long numSamples = frame -> nb_samples;
+  long channels = frame -> channels;
+
+  for (long channel = 0; channel < channels; channel++) {
+    for (long sample = 0; sample < numSamples; sample++) {
+      if (sample + sampleOffset >= startSample
+          and sample + sampleOffset < endSample) {
+
+        double value;
+
+        switch(codecContext -> sample_fmt) {
+          case AV_SAMPLE_FMT_U8:
+            value = ((uint8_t *)(frame -> extended_data))[sample*channels + channel];
+            break;
+          case AV_SAMPLE_FMT_S16:
+            value = ((int16_t *)(frame -> extended_data))[sample*channels + channel];
+            break;
+          case AV_SAMPLE_FMT_S32:
+            value = ((int32_t *)(frame -> extended_data))[sample*channels + channel];
+            break;
+          case AV_SAMPLE_FMT_FLT:
+            value = ((float *)(frame -> extended_data))[sample*channels + channel];
+            break;
+          case AV_SAMPLE_FMT_DBL:
+            value = ((double *)(frame -> extended_data))[sample*channels + channel];
+            break;
+          case AV_SAMPLE_FMT_U8P:
+            value = ((uint8_t *)(frame -> extended_data[channel]))[sample];
+            break;
+          case AV_SAMPLE_FMT_S16P:
+            value = ((int16_t *)(frame -> extended_data[channel]))[sample];
+            break;
+          case AV_SAMPLE_FMT_S32P:
+            value = ((int32_t *)(frame -> extended_data[channel]))[sample];
+            break;
+          case AV_SAMPLE_FMT_FLTP:
+            value = ((float *)(frame -> extended_data[channel]))[sample];
+            break;
+          case AV_SAMPLE_FMT_DBLP:
+            value = ((double *)(frame -> extended_data[channel]))[sample];
+            break;
+          default:
+            throw ProcessingException("Sample format is invalid!");
+            break;
+        }
+
+        output[channel][sample + sampleOffset - startSample] = value;
+      }
+    }
+  }
+
+}
diff --git a/src/AudioSample.cpp b/src/AudioSample.cpp
index 9e6d09f..6ab1cc9 100644
--- a/src/AudioSample.cpp
+++ b/src/AudioSample.cpp
@@ -8,24 +8,43 @@
 #include "SampleSorter/SpectralProcessing.hpp"
 #include "SampleSorter/EqualLoudness.hpp"
 #include "SampleSorter/Octave.hpp"
+#include "SampleSorter/Units.hpp"
 #include "SampleSorter/Tempo.hpp"
 
 AudioSample::AudioSample() {
   tuningCents = 0;
-  tempo = 1;
-  theOne = 0;
 }
 
 AudioSample::AudioSample(
-    std::vector<std::vector<double> > audio, 
-    long sampleRate) {
+    std::vector<std::vector<double> > & audio,
+    long _sampleRate
+    ) : sampleRate(_sampleRate) {
+  totalSeconds = Units::samplesToSeconds(audio[0].size(), sampleRate);
+  tune(audio);
+  findBeat(audio);
+  findChords(audio);
+}
+
+double AudioSample::getTotalSeconds() const {
+  return totalSeconds;
+}
 
-  tune(audio, sampleRate);
-  findBeat(audio, sampleRate);
-  findChords(audio, sampleRate);
+AudioSample::AudioSample(
+    long tuningCents_,
+    double rawBeat,
+    double theOneBin,
+    double totalSeconds_,
+    long sampleRate,
+    std::vector<Octave> chords_) :
+  tempo(rawBeat, theOneBin, sampleRate),
+  sampleRate(sampleRate)
+{
+  totalSeconds = totalSeconds_;
+  tuningCents = tuningCents_;
+  chords = chords_;
 }
 
-void AudioSample::tune(std::vector<std::vector<double> > audio, long sampleRate) {
+void AudioSample::tune(std::vector<std::vector<double> > & audio) {
   std::vector<std::vector<double> > filteredAudio;
   EqualLoudness::filter(filteredAudio, audio, sampleRate);
 
@@ -33,60 +52,93 @@ void AudioSample::tune(std::vector<std::vector<double> > audio, long sampleRate)
   tuningCents = oct.tune();
 }
 
-void AudioSample::findBeat(std::vector<std::vector<double> > audio, long sampleRate) {
-  long hopSize = 1024;
-  long windowRatio = 2;
-
-  std::vector<double> onsets = 
-    SpectralProcessing::onsetEnergy(audio, hopSize, windowRatio);
-  tempo = Tempo::correlationTempo(onsets, hopSize, sampleRate);
-
-  std::pair<double, double> tempoOne = 
-    Tempo::fineTuneTempo(tempo, onsets, 10/60., 0.05/60., hopSize, sampleRate);
-
-  tempo = tempoOne.first;
-  theOne = tempoOne.second;
+void AudioSample::findBeat(std::vector<std::vector<double> > & audio) {
+  double percentageError = 0.5; // 10%
+  double tempoSteps = 1000;
+  double oneSteps = 1000;
+
+  tempo = Tempo(
+      audio, 
+      sampleRate, 
+      percentageError, 
+      tempoSteps,
+      oneSteps
+      );
 }
 
-void AudioSample::findChords(std::vector<std::vector<double> > audio, long sampleRate) {
-  long windowSize = Tempo::tempoToSamples(tempo, sampleRate);
-  std::vector<double> window(windowSize);
-
-  // set up fourier transform
-  long fftSize = windowSize/2 + 1;
-  fftw_complex * fft;
-  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
-  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
-                                           window.data(),
-                                           fft,
-                                           FFTW_ESTIMATE);
-
+void AudioSample::findChords(std::vector<std::vector<double> > & audio) {
+  // filter the audio
   std::vector<std::vector<double> > filteredAudio;
   EqualLoudness::filter(filteredAudio, audio, sampleRate);
   TimeDomainProcessing::unitEnergyPerBeat(
-      filteredAudio, filteredAudio, tempo, sampleRate);
-
-  long theOneSamples = Tempo::secondsToSamples(theOne, sampleRate);
-  long maxWindow = (filteredAudio[0].size() - theOneSamples)/windowSize;
+      filteredAudio, filteredAudio, tempo.getTempo(), sampleRate);
+
+  int windowSize = Units::tempoToSamples(tempo.getTempo(), sampleRate);
+  long startSample = Units::secondsToSamples(tempo.getTheOne(), sampleRate);
+
+  //long k = 0;
+  //double startBeat = 0;
+  //double endBeat = std::min(tempo.getTotalSeconds(), tempo.getKthBeat(k));
+  
+  int numChords = 
+    (Units::secondsToSamples(totalSeconds, sampleRate) - startSample)/windowSize;
+
+  chords.resize(numChords);
+
+  //do {
+    //// set up window
+    //long startSample = Units::secondsToSamples(startBeat, sampleRate);
+    //long endSample = Units::secondsToSamples(endBeat, sampleRate);
+    //long windowSize = endSample - startSample;
+    //if (windowSize > 0 and startSample >= 0) {
+      std::vector<double> window(windowSize);
+
+      // set up fft
+      long fftSize = windowSize/2 + 1;
+      fftw_complex * fft;
+      fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
+      fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windowSize,
+                                               window.data(),
+                                               fft,
+                                               FFTW_ESTIMATE);
+
+      for (int i = 0; i < numChords; i++) {
+
+      Octave chord;
+
+      // for each channel
+      for (long channel = 0; channel < filteredAudio.size(); channel ++) {
+        // window this section
+        //std::cout << "before" << std::endl;
+        for (long i = 0; i < windowSize; i++) {
+          window[i] = filteredAudio[channel][i + startSample];
+          window[i] = window[i] * SpectralProcessing::hammingWindow(i, windowSize);
+        }
+        //std::cout << "after" << std::endl;
+
+        // take fourier transform of this section
+        fftw_execute(fftPlan);
+        // turn it into a chord
+        Octave channelChord(fft, fftSize, 12, sampleRate, tuningCents);
+
+        // add to channel
+        chord.add(chord, channelChord);
+      }
 
-  chords.resize(maxWindow);
 
-  for (long hop = 0; hop < maxWindow; hop++) {
-    for (long channel = 0; channel < filteredAudio.size(); channel ++) {
-      for (long i = 0; i < windowSize; i++) {
-        window[i] = filteredAudio[channel][theOneSamples + i + hop * windowSize];
-        window[i] = window[i] * SpectralProcessing::hammingWindow(i, windowSize);
+      // add to chords
+      chords[i] = chord;
       }
 
-      fftw_execute(fftPlan);
-      Octave channelOctave(fft, fftSize, 12, sampleRate, tuningCents);
 
-      chords[hop].add(chords[hop], channelOctave);
-    }
-  }
+      fftw_free(fft);
+      fftw_destroy_plan(fftPlan);
+    //}
 
-  fftw_free(fft);
-  fftw_destroy_plan(fftPlan);
+    //k += 1;
+    //startBeat = endBeat;
+    //endBeat = std::min(totalSeconds, tempo.getKthBeat(k));
+  //} while (startBeat != tempo.getTotalSeconds());
 }
 
 long AudioSample::getTuningCents() const {
@@ -96,17 +148,20 @@ double AudioSample::getTuningCentsFreqRatio() const {
   return std::pow(2., tuningCents/1200.);
 }
 double AudioSample::getBeatRaw() const {
-  return tempo;
+  return tempo.getTempo();
 }
 double AudioSample::getBeatWithTuning() const {
-  return tempo * getTuningCentsFreqRatio();
+  return tempo.getTempo() * getTuningCentsFreqRatio();
 }
 double AudioSample::getTheOneRaw() const {
-  return theOne;
+  return tempo.getTheOne();
 }
 double AudioSample::getTheOneWithTuning() const {
-  return theOne/getTuningCentsFreqRatio();
+  return tempo.getTheOne()/getTuningCentsFreqRatio();
 }
 std::vector<Octave> AudioSample::getChords() const {
   return chords;
 }
+long AudioSample::getSampleRate() const {
+  return sampleRate;
+}
diff --git a/src/Octave.cpp b/src/Octave.cpp
index d0be445..0daed0f 100644
--- a/src/Octave.cpp
+++ b/src/Octave.cpp
@@ -110,11 +110,11 @@ void Octave::addPeak(double peakFreq, double peakValue, double baseOffset) {
     double rightAmp = rightPercent * peakValue;
     double leftAmp = leftPercent * peakValue;
 
-    if (std::isnan(rightAmp) or std::isnan(leftAmp) ) {
-      std::cout << "ahhhh" << std::endl;
-      std::cout << peakValue << std::endl;
-      std::cout << desiredBin << std::endl;
-    }
+    //if (std::isnan(rightAmp) or std::isnan(leftAmp) ) {
+      //std::cout << "ahhhh" << std::endl;
+      //std::cout << peakValue << std::endl;
+      //std::cout << desiredBin << std::endl;
+    //}
 
     spectrogram[leftBin] += leftAmp;
     spectrogram[rightBin] += rightAmp;
diff --git a/src/Plotting.cpp b/src/Plotting.cpp
index d7692f4..01391bc 100644
--- a/src/Plotting.cpp
+++ b/src/Plotting.cpp
@@ -20,3 +20,21 @@ void Plotting::plotPair(std::vector<std::pair<double, double> > xy) {
   gp << "plot" << gp.file1d(xy) << "w l" << std::endl;
   std::cin.get();
 }
+
+void Plotting::plotLineAndMarkers(
+    std::vector<std::pair<double, double> > line,
+    std::vector<double> markers,
+    double pointHeight) {
+
+  std::vector<std::pair<double, double> > markerPoints(markers.size());
+  for (long i = 0; i < markerPoints.size(); i++) {
+    markerPoints[i] = std::make_pair(markers[i], pointHeight);
+  }
+
+  Gnuplot gp;
+  gp << "plot '-' w l, '-' w p" << std::endl;
+  gp.send1d(line);
+  gp.send1d(markerPoints);
+  std::cin.get();
+}
+
diff --git a/src/PythonExternals.cpp b/src/PythonExternals.cpp
index 8f3fbdb..ffd4ebb 100644
--- a/src/PythonExternals.cpp
+++ b/src/PythonExternals.cpp
@@ -9,9 +9,14 @@
 extern "C" {
   SampleFile * NewAbletonSampleFile(const char * filePath) {
     AbletonSampleFile * s = new AbletonSampleFile(filePath);
-    s -> process();
     return s;
   }
+  bool process(SampleFile * s) {
+    return s -> process();
+  }
+  const char * getFileName(SampleFile * s) {
+    return s -> getFileName().c_str();
+  }
   long getTuningCents(SampleFile * s) {
     return s -> getAudioSample() -> getTuningCents();
   }
@@ -35,6 +40,9 @@ extern "C" {
     }
     return chordsArray;
   }
+  void writeToFile(AbletonSampleFile * s) {
+    s -> writeToFile();
+  }
   void deleteChords(double ** chords, size_t chordsSize) {
     for (size_t i = 0; i < chordsSize; i++) {
       delete chords[i];
diff --git a/src/SampleFile.cpp b/src/SampleFile.cpp
index ab8ce71..2e4b3eb 100644
--- a/src/SampleFile.cpp
+++ b/src/SampleFile.cpp
@@ -5,13 +5,24 @@
 
 #include "SampleSorter/SampleFile.hpp"
 #include "SampleSorter/AudioSample.hpp"
+#include "SampleSorter/ProcessingException.hpp"
 
 SampleFile::SampleFile(std::string filePath_) {
   filePath = filePath_;
 }
 
-void SampleFile::process() {
-  sample = AudioSample(getWaves(), getSampleRate());
+bool SampleFile::process() {
+  try {
+    if (not readMetaData()) {
+      long sampleRate;
+      std::vector< std::vector<double> > waves = extractAudio(&sampleRate);
+      sample = AudioSample(waves, sampleRate);
+    }
+    return true;
+  } catch (ProcessingException & e) {
+    std::cout << e.getMessage() << std::endl;
+    return false;
+  }
 }
 
 std::string SampleFile::getFilePath() {
@@ -25,55 +36,3 @@ std::string SampleFile::getFileName() {
 AudioSample * SampleFile::getAudioSample() {
   return &sample;
 }
-
-//void SampleFile::writeToMIDI(std::string fileName) {
-  //// Find max
-  //double max = 0;
-  //for (long i = 0; i < chords.size(); i++) {
-    //for (long bin = 0; bin < 12; bin ++) {
-      //getSpectrogram;
-      //max = std::max(max, spec[bin]);
-    //}
-  //}
-
-  //ratio = 127./max;
-
-  //// open midi file
-  //// set tempo
-  //// add track
-  //// for every chord
-  //// for every bin
-  //// add an event
-//}
-
-
-
-
-//bool Sample::isCompatible(const Sample & other) {
-  //// get larger tempo
-  //double larger, smaller;
-  //larger = std::max(getBeatWithTuning(), other.getBeatWithTuning());
-  //smaller = std::min(getBeatWithTuning(), other.getBeatWithTuning());
-  //long ratio = std::round(larger/smaller);
-  //double tuningSteps = 12 * std::log2(larger/(ratio * smaller));
-  ////double integerDeviation = std::fmod(tuningSteps, 1.);
-  //double integerDeviation = tuningSteps - std::round(tuningSteps);
-
-  //if (std::abs(integerDeviation) < 0.05) {
-    //std::cout << "'" << getName() << "' is compatible with '" << other.getName() <<"'" << std::endl;
-    //std::cout << "Tempos: " << 60*getBeatWithTuning() <<", " << 60*other.getBeatWithTuning() << std::endl;
-    //std::cout << "Tuning cents: " << tuningCents <<", " << other.tuningCents << std::endl;
-    //std::cout << "Tuning ratio: " << std::round(tuningSteps) << " steps and " << integerDeviation * 100 << " cents" << std::endl << std::endl;
-  //}
-
-  //return integerDeviation < 0.05;
-//}
-
-  // find whether that ratio is an even number of cents
-
-//bool Sample::isSimilar(const & Sample other) {
-  //// coarse tune so their tempos are equal
-  //// or are integer ratios of each other
-  //// get similarity of all octaves.
-  //// Look for large chains
-//}
diff --git a/src/Tempo.cpp b/src/Tempo.cpp
index ef8bdd8..70dde29 100644
--- a/src/Tempo.cpp
+++ b/src/Tempo.cpp
@@ -1,168 +1,260 @@
 #include <vector>
-#include <iostream>
+#include <cmath>
+#include <random>
+
+#include <fftw3.h>
 
-#include "Plotting/Plotting.hpp"
 #include "SampleSorter/Tempo.hpp"
+#include "SampleSorter/Units.hpp"
 #include "SampleSorter/SpectralProcessing.hpp"
+#include "Plotting/Plotting.hpp"
 
-double Tempo::tempoToSeconds(double tempo) {
-  return 1/tempo;
-}
+const long Tempo::HOP_SIZE = 1024;
 
-double Tempo::tempoToSamples(double tempo, long sampleRate) {
-  return secondsToSamples(tempoToSeconds(tempo), sampleRate);
+Tempo::Tempo() : sampleRate(1) {
+  tempo = 1;
+  theOneBin = 0;
 }
 
-double Tempo::secondsToSamples(double seconds, long sampleRate) {
-  return seconds * sampleRate;
+double Tempo::getTheOneBin() const {
+  return theOneBin;
 }
 
-double Tempo::samplesToSeconds(double samples, long sampleRate) {
-  return samples/double(sampleRate);
+double Tempo::getTheOne() const {
+  return Units::binsToSeconds(theOneBin, HOP_SIZE, sampleRate);
 }
 
-double Tempo::tempoToBins(double tempo, long hopSize, long sampleRate) {
-  return sampleRate/(tempo * hopSize);
+Tempo::Tempo(double tempo_, double theOneBin_, long sampleRate_) : sampleRate(sampleRate_) {
+  tempo = tempo_;
+  theOneBin = theOneBin_;
 }
 
-double Tempo::binsToSamples(long bin, long hopSize) {
-  return bin * hopSize;
-}
+void Tempo::fineTuneTempo(
+    const double percentageError,
+    const int steps,
+    const std::vector<double> & onsets
+    ) {
+  double guessTempo = tempo;
+
+  double minValue = onsets.size();
+  double minTempo = guessTempo;
+  int minOneBin = 0;
+
+  // Guess a tempo in the range
+  // iterate over bins, because
+  // the error will be linear with number of bins
+  double midBins = Units::tempoToBins(guessTempo, HOP_SIZE, sampleRate);
+  double minBins = midBins * (1 - percentageError);
+  double maxBins = midBins * (1 + percentageError);
+  for (double bins = minBins;
+      bins < maxBins;
+      bins += (maxBins - minBins)/double(steps)) {
+
+    tempo = Units::binsToTempo(bins, HOP_SIZE, sampleRate);
+    
+    // guess a one in the range
+    for (theOneBin = 0;
+        theOneBin < bins;
+        theOneBin ++) {
+      // calculate the value
+      double newValue = getValue(onsets, true);
+
+      // choose one that minimize energy off of beat
+      if (newValue < minValue) {
+        minValue = newValue;
+        minTempo = tempo;
+        minOneBin = theOneBin;
+      }
+    }
+
+  }
 
-double Tempo::binsToSeconds(long bin, long hopSize, long sampleRate) {
-  return samplesToSeconds(binsToSamples(bin, hopSize), sampleRate);
+  tempo = minTempo;
+  theOneBin = minOneBin;
 }
 
+void Tempo::fineTuneTheOne(
+    const std::vector<double> & onsets,
+    const int steps) {
+  // the first beat
+  double beatOneBins = Units::tempoToBins(tempo, HOP_SIZE, sampleRate);
 
-double Tempo::tempoValue(double tempo, 
-                         double theOne,
-                         std::vector<double> onsets, 
-                         long hopSize,
-                         long sampleRate,
-                         bool bidirectional) {
-  double value = 0;
-  for (long i = 0; i < onsets.size(); i++) {
-    double distanceFromBeat = (i - theOne)/tempoToBins(tempo, hopSize, sampleRate);
-    if (bidirectional) {
-      distanceFromBeat = std::abs(distanceFromBeat - std::round(distanceFromBeat));
-    } else {
-      distanceFromBeat = distanceFromBeat - std::floor(distanceFromBeat);
+  double minValue = onsets.size();
+  double minOneBin = 0;
+
+  // for all possible ones
+  for (theOneBin = 0;
+      theOneBin < beatOneBins;
+      theOneBin += beatOneBins/double(steps)) {
+    double value = getValue(onsets, false);
+
+    if (value < minValue) {
+      minValue = value;
+      minOneBin = theOneBin;
     }
+  }
+  theOneBin = minOneBin;
+}
+
+double Tempo::getTempo() const {
+  return tempo;
+}
 
-    value += distanceFromBeat * onsets[i];
+double Tempo::distanceFromBeat(
+    double bin,
+    bool bidirectional,
+    long hopSize,
+    long sampleRate
+    ) const {
+  double distance = 
+    (bin - theOneBin)/Units::tempoToBins(tempo, hopSize, sampleRate);
+  if (bidirectional) {
+    distance = std::abs(distance - std::round(distance));
+  } else {
+    distance = distance - std::floor(distance);
   }
+  return distance;
+}
+
+double Tempo::getValue(
+    const std::vector<double> & onsets,
+    bool bidirectional
+    ) const {
 
-  return value;//tempoToBins(tempo, hopSize, sampleRate);
+  double value = 0;
+  for (size_t bin = 0; bin < onsets.size(); bin ++) {
+    value += distanceFromBeat(bin, HOP_SIZE, sampleRate, bidirectional) * onsets[bin];
+  }
+  
+  return value;
 }
 
+Tempo::Tempo(
+    const std::vector< std::vector<double> > & audio,
+    long sampleRate_,
+    double percentageError,
+    int tempoSteps,
+    int oneSteps) : sampleRate(sampleRate_) {
+
+  std::cout << audio.size() << ", " << audio[0].size() << std::endl;
+  std::cout << sampleRate << std::endl;
 
-double Tempo::correlationTempo(std::vector<double> onsets, 
-                               long hopSize, 
-                               long sampleRate) {
-  //Plotting::plotVector(onsets, hopSize);
+  //std::cout << "starting tempo" << std::endl;
+
+  long windowRatio = 2;
+  std::vector<double> onsets =
+    SpectralProcessing::onsetEnergy(audio, HOP_SIZE, windowRatio);
+
+  findCorrelationTempo(onsets);
+
+  std::cout << "guessed tempo: " << getTempo() * 60. << std::endl;
+
+  fineTuneTempo(
+      percentageError,
+      tempoSteps,
+      onsets
+      );
+
+  std::cout << "fine tuned tempo: " << getTempo() * 60. << std::endl;
+
+  //aCoefficients.resize(degrees);
+  //bCoefficients.resize(degrees);
+  //std::fill(aCoefficients.begin(), aCoefficients.end(), 0);
+  //std::fill(bCoefficients.begin(), bCoefficients.end(), 0);
+  //gradientDescent(onsets, hopSize, sampleRate);
+
+  std::cout << "the One: " << getTheOne() << std::endl;
+
+  fineTuneTheOne(onsets, oneSteps);
+
+  std::cout << "the fine tuned One: " << getTheOne() << std::endl;
+
+  //for (long i = 0; i < getDegree(); i++) {
+    //std::cout << "a[" << i << "]=" << aCoefficients[i] << std::endl;
+    //std::cout << "a[" << i << "]=" << bCoefficients[i] << std::endl;
+  //}
+
+  plotOnsetsWithBeats(onsets);
+  //plotBeats(0.001);
+  //plotTempo(0.001);
+  //
+
+}
+
+void Tempo::findCorrelationTempo(
+    const std::vector<double> & onsets
+    ) {
 
   // filter out frequencies less than 2 per clip
   // Clips with tempo must contain at least 2 beats
-  double highPass = 2./binsToSeconds(onsets.size(), hopSize, sampleRate);
+  double highPass = 2./Units::binsToSeconds(onsets.size(), HOP_SIZE, sampleRate);
   // Also beat cannot be below 1/8 persecond
   highPass = std::max(highPass, 1.);
 
   // No tempos will exist above 1000bpm
   double lowPass = 1000/60.;
 
-  std::vector<double> correlation = SpectralProcessing::autoCorrelation(onsets,
-                                                    sampleRate/hopSize,
-                                                    highPass,
-                                                    lowPass);
+  std::vector<double> correlation = 
+    SpectralProcessing::autoCorrelation(
+        onsets,
+        sampleRate/HOP_SIZE,
+        highPass,
+        lowPass);
 
   // only take first half of correlation
   long fftSize = correlation.size()/4 + 1;
   fftw_complex * fft;
   fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
-  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(correlation.size()/2,
-                                              correlation.data(),
-                                              fft,
-                                              FFTW_ESTIMATE);
+  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(
+      correlation.size()/2,
+      correlation.data(),
+      fft,
+      FFTW_ESTIMATE);
 
   //Window the correlation
   for (long i = 0; i < correlation.size()/2; i++) {
-    correlation[i] = correlation[i] * SpectralProcessing::hammingWindow(i, correlation.size()/2);
+    correlation[i] = correlation[i] 
+              * SpectralProcessing::hammingWindow(i, correlation.size()/2);
   }
 
   fftw_execute(fftPlan);
   fftw_destroy_plan(fftPlan);
 
   std::vector<std::pair<double, double> > peaks = 
-    SpectralProcessing::findPeaks(fft, fftSize, sampleRate/hopSize);
+    SpectralProcessing::findPeaks(fft, fftSize, sampleRate/HOP_SIZE);
 
-  //Plotting::plotPair(peaks);
+  Plotting::plotPair(peaks);
 
   fftw_free(fft);
 
   double max = 0;
-  double tempo = 0;
+  tempo = 0;
   for (long i = 0; i < peaks.size(); i++) {
     if (peaks[i].second > max) {
       max = peaks[i].second;
       tempo = peaks[i].first;
     }
   }
-
-  return tempo;
 }
 
-std::pair<double, double> Tempo::fineTuneTempo(double tempo, 
-    std::vector<double> onsets, 
-    double BPSrange, 
-    double BPSstepSize,
-    long hopSize,
-    long sampleRate) {
-
-  std::vector<double> tempoValues;
-
-  //Plotting::plotVector(onsets);
-
-  // first find the bidirectional minima
-  double minTempoValue = onsets.size();
-  double minTempo = tempo;
-  for (double tempoAdjustment = -BPSrange; 
-       tempoAdjustment <=BPSrange; 
-       tempoAdjustment += BPSstepSize) {
-    double trialTempo = tempo + tempoAdjustment;
-    // adjust it by size of beats to bins -> less bins -> easier
-    double thisMinValue = onsets.size();
-    double thisMinTempo = 0;
-    for (double theOne = 0; 
-         theOne < tempoToBins(trialTempo, hopSize, sampleRate); 
-         theOne ++) {
-      double newTempoValue = tempoValue(trialTempo, theOne, onsets, hopSize, sampleRate, true);
-      if (newTempoValue < minTempoValue) {
-        minTempoValue = newTempoValue;
-        minTempo = trialTempo;
-      }
+// plot onsets with beats
+void Tempo::plotOnsetsWithBeats(
+    const std::vector<double> & onsets
+    ) const {
 
-      if (newTempoValue < thisMinValue) {
-        thisMinValue = newTempoValue;
-        thisMinTempo = trialTempo;
-      }
-    }
-    tempoValues.push_back(thisMinValue);
+  std::vector<double> beats;
+  double beat = getTheOne();
+  while (beat < Units::binsToSeconds(onsets.size(), HOP_SIZE, sampleRate)) {
+    beats.push_back(beat);
+    beat += Units::tempoToSeconds(tempo);
   }
 
-  // then search for the start of the beat
-  minTempoValue = onsets.size();
-  std::vector<double> minValues;
-  double minOne = 0;
-  for (double theOne = 0; 
-       theOne < tempoToBins(minTempo, hopSize, sampleRate)+1; 
-       theOne ++) {
-    double newTempoValue = tempoValue(minTempo, theOne, onsets, hopSize, sampleRate, false);
-    if (newTempoValue < minTempoValue) {
-      minTempoValue = newTempoValue;
-      minOne = theOne;
-    }
-    minValues.push_back(newTempoValue);
+  std::vector<std::pair<double, double> > onsetSeconds(onsets.size());
+
+  for (long bin = 0; bin < onsets.size(); bin++) {
+    onsetSeconds[bin] = std::make_pair(Units::binsToSeconds(bin, HOP_SIZE, sampleRate), onsets[bin]);
   }
 
-  return std::make_pair(minTempo, binsToSeconds(minOne, hopSize, sampleRate));
+  Plotting::plotLineAndMarkers(onsetSeconds, beats, 0.5);
 }
+
diff --git a/src/Units.cpp b/src/Units.cpp
index 8c6da48..d4a8dfe 100644
--- a/src/Units.cpp
+++ b/src/Units.cpp
@@ -15,6 +15,10 @@ double Units::secondsToSamples(double seconds, long sampleRate) {
   return seconds * sampleRate;
 }
 
+double Units::secondsToTempo(double seconds) {
+  return 1/seconds;
+}
+
 double Units::samplesToSeconds(double samples, long sampleRate) {
   return samples/double(sampleRate);
 }
@@ -23,10 +27,22 @@ double Units::tempoToBins(double tempo, long hopSize, long sampleRate) {
   return sampleRate/(tempo * hopSize);
 }
 
-double Units::binsToSamples(long bin, long hopSize) {
+double Units::binsToSamples(double bin, long hopSize) {
   return bin * hopSize;
 }
 
-double Units::binsToSeconds(long bin, long hopSize, long sampleRate) {
+double Units::binsToSeconds(double bin, long hopSize, long sampleRate) {
   return samplesToSeconds(binsToSamples(bin, hopSize), sampleRate);
 }
+
+double Units::binsToTempo(double bin, long hopSize, long sampleRate) {
+  return secondsToTempo(binsToSeconds(bin, hopSize, sampleRate));
+}
+
+double Units::samplesToBins(double samples, long hopSize) {
+  return samples/double(hopSize);
+}
+
+double Units::secondsToBins(double seconds, long hopSize, long sampleRate) {
+  return samplesToBins(secondsToSamples(seconds, sampleRate), hopSize);
+}
diff --git a/test/MPEG_test.cpp b/test/MPEG_test.cpp
new file mode 100644
index 0000000..d7910b2
--- /dev/null
+++ b/test/MPEG_test.cpp
@@ -0,0 +1,26 @@
+#include <string>
+#include <vector>
+
+#include <SampleSorter/AudioFile.hpp>
+#include <Plotting/Plotting.hpp>
+
+#include "gtest/gtest.h"
+
+TEST(AudioFileTest, BasicRead) {
+  std::string filePath = "../TestFiles/MPEGFiles/Aintnocrime.mp3";
+  //std::string filePath = "../TestFiles/MPEGFiles/ComeTogether.m4a";
+  //std::string filePath = "../TestFiles/MPEGFiles/WouldntItBeNice.mp3";
+  //std::string filePath = "../TestFiles/MPEGFiles/ArcadeFire.mp3";
+  //std::string filePath = "../TestFiles/MPEGFiles/Bear.mp3";
+  //
+  long sampleRate;
+  
+  std::vector<std::vector<double> > output = AudioFile::read(filePath, 0, 10, &sampleRate);
+
+  Plotting::plotVector(output[0]);
+}
+
+int main(int argc, char ** argv) {
+  ::testing::InitGoogleTest(&argc, argv);
+  return RUN_ALL_TESTS();
+}
diff --git a/test/tensorFlow_test.cpp b/test/tensorFlow_test.cpp
new file mode 100644
index 0000000..eeb626d
--- /dev/null
+++ b/test/tensorFlow_test.cpp
@@ -0,0 +1,7 @@
+#include <iostream>
+
+//#include <tensorflow>
+
+int main() {
+  std::cout << "yeah, I do work" << std::endl;
+}
