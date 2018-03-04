#include <string>
#include <iostream>

#include <SampleSorter/AbletonSampleFile.hpp> 
#include <SampleSorter/EqualLoudness.hpp> 
#include <SampleSorter/SpectralProcessing.hpp>
#include <SampleSorter/Units.hpp>
#include <Plotting/Plotting.hpp>

std::vector<double> getFFTMag(std::vector<double> audio) {
  long fftSize = audio.size()/2 + 1;

  std::vector<double> windows(audio.size());
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(windows.size(),
                                 windows.data(),
                                 fft,
                                 FFTW_ESTIMATE |
                                 FFTW_DESTROY_INPUT);

  // Window the audio
  SpectralProcessing::hammingWindow(windows, audio);

  // Take the fourier transform
  fftw_execute(fftPlan);

  std::vector<double> output;
  output.resize(fftSize);
  for (long i = 0; i < fftSize; i++) {
     output[i] = std::abs(SpectralProcessing::fftToComplex(fft, i));
  }

  return output;
}

int main(int argc, char ** argv) {
  // The test sample file
  std::string sampleFile = "../test/User Library/Samples/sample.alc";
  // The user library root
  std::string userLibrary = "../test/User Library/";

  // Open the sample file and process it to get file name, etc.
  AbletonSampleFile sample(sampleFile, userLibrary, true);
  sample.process();

  // Get the sample rate and audio
  long sampleRate;
  std::vector< std::vector<double> > audio = sample.extractAudio(&sampleRate);
  // Plot the audio
  double samplesToSeconds = Units::samplesToSeconds(1, sampleRate);
  std::cout << "Plotting audio." << std::endl;
  Plotting::plotVector(audio[0], "Audio", "Time (seconds)", "Amplitude", samplesToSeconds);
  std::cout << "Plotting audio FFT." << std::endl;
  std::vector<double> audioFFTMag = getFFTMag(audio[0]);
  double binToFreq = sampleRate/double(audio[0].size());
  Plotting::plotVector(audioFFTMag, "Audio FFT", "Frequency (hz)", "Amplitude", binToFreq, 0, false, 10);
  
  // Filter by equal loudness
  std::vector<std::vector<double> > filteredAudio;
  EqualLoudness::filter(filteredAudio, audio, sampleRate);
  std::cout << "Plotting filtered audio." << std::endl;
  // Plot the filtered audio
  Plotting::plotVector(filteredAudio[0], "A-Weighted Audio", "Time (seconds)", "Amplitude", samplesToSeconds);
  std::cout << "Plotting filtered audio FFT." << std::endl;
  std::vector<double> filteredAudioFFTMag = getFFTMag(filteredAudio[0]);
  Plotting::plotVector(filteredAudioFFTMag, "A-Weighted Audio FFT", "Time (seconds)", "Amplitude", binToFreq, 0, false, 10);
  
  // Convert to octave
  Octave oct(filteredAudio, 1200., sampleRate);
  // Plot it
  std::cout << "Plotting octave." << std::endl;
  oct.plot("Octave-Wrapped Spectrum", "12-Tone Equal Tempered Pitch (0 = A)", "Amplitude");

  // Tune the octave
  long tuningCents = oct.tune();
  std::cout << "Tuned by " << tuningCents << " cents." << std::endl;
  // Plot it
  oct.plot("Tuned Octave-Wrapped Spectrum", "12-Tone Equal Tempered Pitch (0 = A)", "Amplitude");

  // Convert to a discretized octave
  Octave octSmall(filteredAudio, 12., sampleRate, tuningCents = tuningCents);
  // Plot it
  std::cout << "Plotting octave." << std::endl;
  octSmall.plot("Discretized Octave-Wrapped Spectrum", "12-Tone Equal Tempered Pitch (0 = A)", "Amplitude", true);

  // Get the onset energy
  long hopSize = 1024;
  long windowRatio = 2;
  std::vector<double> onsets =
    SpectralProcessing::onsetEnergy(audio, hopSize, windowRatio);
  // Plot
  double binsToSeconds = Units::binsToSeconds(1, hopSize, sampleRate);
  std::cout << "Plotting onset energy" << std::endl;
  Plotting::plotVector(onsets, "Onset Energy", "Time (seconds)", "Amplitude", binsToSeconds);
  
  // Find the auto correlation
  double highPass = 2./Units::binsToSeconds(onsets.size(), hopSize, sampleRate);
  highPass = std::max(highPass, 1.);
  std::vector<double> correlation = 
    SpectralProcessing::autoCorrelation(
        onsets,
        sampleRate/hopSize,
        highPass,
        1000/60);
  std::cout << "Plotting autocorrelation." << std::endl;
  Plotting::plotVector(correlation, "Onset Energy Autocorrelation", "Time (seconds)", "Amplitude", binsToSeconds);
  
  // only take first half of correlation
  long fftSize = correlation.size()/4 + 1;
  fftw_complex * fft;
  fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
  fftw_plan fftPlan = fftw_plan_dft_r2c_1d(
      correlation.size()/2,
      correlation.data(),
      fft,
      FFTW_ESTIMATE);

  //Window the correlation
  for (long i = 0; i < correlation.size()/2; i++) {
    correlation[i] = correlation[i] 
              * SpectralProcessing::hammingWindow(i, correlation.size()/2);
  }

  fftw_execute(fftPlan);
  fftw_destroy_plan(fftPlan);
  std::vector<double> correlationFFT;
  correlationFFT.resize(fftSize);
  for (long i = 0; i < fftSize; i++) {
     correlationFFT[i] = std::abs(SpectralProcessing::fftToComplex(fft, i));
  }
  double binsToBPM = 60. * (sampleRate/double(hopSize))/(2*fftSize + 1);
  std::cout << "Plotting autocorrelation FFT." << std::endl;
  Plotting::plotVector(correlationFFT, "Onset Energy Autocorrelation FFT", "Tempo (BPM)", "Amplitude", binsToBPM, 0, false, 50);

  // Fine tune
  Tempo tempo = Tempo(
      audio, 
      sampleRate, 
      0.2,
      2000,
      1000
      );
  std::cout << "Tempo: " << 60. * tempo.getTempo() << std::endl;
  std::cout << "Plotting onset energy with tempo." << std::endl;
  // Plot
  tempo.plotOnsetsWithBeats(onsets, "Onset Energy with Tuned Tempo", "Time (seconds)", "Amplitude");
  std::cout << "Tempo with tuning: " << 60. * (sample.getAudioSample() -> getBeatWithTuning()) << std::endl;
}
