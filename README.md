# SampleSorter

This code analyzes a library of audio samples and determines each sample's tempo, tuning, and harmonic profile.
Using these values we determine pairs of samples that, when played at particular speeds, 
have the same tempo, have the same tuning and are harmonically coherent.
It is intended as a compositional tool for making sample-based music.

## Motivation

The [sampling](https://en.wikipedia.org/wiki/Sampling_(music)) technique this work serves to aid is sometimes known as "[plunderphonics](https://en.wikipedia.org/wiki/Plunderphonics)" or a "sound collage".
It was popularized by The Avalanches in their masterpiece [Since I Left You](https://www.youtube.com/watch?v=LhBacKKEyBU).
This technique involves layering many samples on top of each other.
DJs employ a related technique called [harmonic mixing](https://en.wikipedia.org/wiki/Harmonic_mixing) to transition smoothly between different songs.
Even though the samples are often untouched (no chopping or added effects), 
the combination of different moods and genres and in which they are placed creates an entirely new sound.
There is something incredibly moving and surreal in being able to take recordings from different genres, different decades, and different continents and realizing them as different parts of the same song.

Despite its beauty, music made with this technique is rare, perhaps because of how difficult it is to create;
the process of finding groups of samples that can be played simultaneously is extremely time consuming.
This is partly because two fundamental musical properties, tempo and pitch, are inherently linked.
Lowering or raising the pitch of a song slows down or speeds up the tempo.
If two samples do in fact have the same tempo and tuning 
(tuning measured in cents as the deviation from [12-tone equal tempermant](https://en.wikipedia.org/wiki/Equal_temperament) with [A4 = 440hz](https://en.wikipedia.org/wiki/A440_(pitch_standard))),
they might not be harmonically coherent.
The samples could be in conflicting keys or modes, again greatly reducing the number valid sample pairs.
Once a musician has finally found samples that have consistent rhythm, tuning and harmony,
they then get to choose which samples should actually go together based on aesthetics and composition.

It should be noted that a process called [time stretching](https://en.wikipedia.org/wiki/Audio_time_stretching_and_pitch_scaling) is able to change the speed of audio without affecting the pitch.
However it necessitates windowing the audio which introduces digital artifacts particularly at transients.
Modern algorithms cover up most of these artifacts, but there is nevertheless something uncanny about time-stretched audio even for small changes in speed.
Even more advanced algorithm's like the technology that exists in the commericial software [Melodyne](https://en.wikipedia.org/wiki/Celemony_Software#DNA_Direct_Note_Access) allow for the pitches of individual notes within a single audio sample to be modified independenly; changing the internal harmony. But once again this can cause harsh digital artifacts.

## Solution Overview

This code works to make the sampling process easier by automatically detecting pairs of samples that can be played at the same time and still sound "musical".

I have a collection of audio samples in Ableton's [\*.alc file format](https://help.ableton.com/hc/en-us/articles/209769625-Live-specific-file-types-adg-als-alp-) (a convenient way to make references to sections of audio).
These are small (~3-30s long) snippets of music from a variety of genres; 
soul, jazz, funk, rock, hip-hop, disco, electronic, experimental, etc.
My own collection is not public due to copyright, 
but there is a large online database of the samples used in popular music at [WhoSampled.com](https://www.whosampled.com/most-sampled-tracks/)

Each of these samples is analyzed using spectral analysis to determine it's tempo, tuning and harmonic profile. The harmonic profile encodes how prevalent each of the [12 notes](https://en.wikipedia.org/wiki/Equal_temperament) is in the sample. Since there is a closed-form relation between tempo and pitch, we can easily determine which pairs of samples can be repitched to have the same tempo and still be in tune. A song with tempo 
![tempo_1](https://latex.codecogs.com/gif.latex?\tau_1)
needs to be retuned by 
![tuning](https://latex.codecogs.com/gif.latex?\delta_{12})
[cents](https://en.wikipedia.org/wiki/Cent_(music)) in order to have tempo
![tempo_1](https://latex.codecogs.com/gif.latex?\tau_2):

![tuning_equation](https://latex.codecogs.com/gif.latex?\delta_{12}=1200\cdot\log_2\left(\frac{\tau_2}{\tau_1}\right))

So long as the tuning is approximately an integer number of [semitones](https://en.wikipedia.org/wiki/Semitone),
![tempo_1](https://latex.codecogs.com/gif.latex?\left|\delta_{12}\right|\mod{50}<\epsilon),
then the pair is "in tune."

The harmonic profiles of these tuned sample pairs are combined and fed into a classifier which rates them based on how harmonically coherent they are. For example a profile which is dominated by the notes "C" "E" and "G" (a major chord) would be labeled "coherent", where as a profile where all the notes have equal weight would be labeled "incoherent". The classifier, which is a fully connected neural network, is trained on the sample collection itself using positive unlabeled learning.

## Solution Description

![audio](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/audio.png)

### Tempo

To detect the tempo we first compute the 
[onsets](https://en.wikipedia.org/wiki/Onset_(audio))
of the audio using the complex-domain method described in
[On the Use of Phase and Energy for Musical Onset Detection in the Complex Domain](https://www.researchgate.net/profile/Mark_Sandler2/publication/3343056_On_the_Use_of_Phase_and_Energy_for_Musical_Onset_Detection_in_the_Complex_Domain/links/5412b6110cf2bb7347dafd25/On-the-Use-of-Phase-and-Energy-for-Musical-Onset-Detection-in-the-Complex-Domain.pdf).
The onset signal is high when there are changes in amplitude or frequency, which typically represent beats:

![onset_energy](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/onset_energy.png)

We then take the [autocorrelation](https://en.wikipedia.org/wiki/Autocorrelation) of the onset signal which reveals the periodic signal in the onsets:
![onset_energy_autocorrelation](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/onset_energy_autocorrelation.png)

The fourier transform of the autocorrelation is peaked around the tempo! We find the peak using the algorithm described in [Estimation of Frequency, Amplitude, and Phase from the DFT of a Time Series](https://pdfs.semanticscholar.org/df2e/2b3ae9d784e19ea0840f8bb26ff622b17c22.pdf):
![onset_energy_autocorrelation_fft](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/onset_energy_autocorrelation_fft.png)

We then fine tune the peak tempo by iterating over a small window around that value.
In each iteration we perform a weighted sum of the onsets. The weights are low close to each beat, as determined by the tempo in that particular iteration (marked in blue below), and high away from them.
By choosing the tempo that minimizes that sum, the tempo very clearly aligns with the onsets:
![onset_energy_with_tempo](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/onset_energy_with_tempo.png)

### Octave-Mapped Spectrogram

![audio_fft](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/audio_fft.png)
![audio_filtered_fft](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/audio_filtered_fft.png)
![octave](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/octave.png)
![octave_tuned](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/octave_tuned.png)
![octave_discretized](https://raw.githubusercontent.com/sportdeath/SampleSorter/master/media/octave_discretized.png)

### Classification

Using the octave-mapped spectrogram as a compact profile of a sample's harmonic content,
we can now classify pairs of samples as being harmonic or inharmonic.
While throughout history many "[rules of harmony](https://en.wikipedia.org/wiki/Harmony#Historical_rules)" have been proposed, those rules are often broken and change over time. Moreover they would be extremely difficult to impliment consistently.
Instead we take advantage of the fact that we already have a dataset of samples which are "harmonic".

Using this dataset, classifying pairs of samples is a positive-unlabeled learning problem.
All of the samples in the dataset are considered to be harmonic examples.
We can combine any pair of these samples to create an unlabeled example.

I have chosen to classify these samples using a multilayer fully connected neural network.
Since the network has so many free parameters, PU-learning is prone to overfit so we use a non-negative loss function:

[Latexed function]

In addition we normalize all the inputs, add L2-regularization, batch normalization, and dropout to prevent overfitting. We also augment the data by rotating the octaves uniformly at random. To have better rotational invariance we employ a variant of the [spatial transoformer netowork](link).

[GIF of octave transform]

## Results

[Examples of found pairs]

## Dependencies

- Cmake
- ffmpeg
- tinyxml2
- fftw3
- boost

## References

- [On the Use of Phase and Energy for Musical Onset Detection in the Complex Domain](https://www.researchgate.net/profile/Mark_Sandler2/publication/3343056_On_the_Use_of_Phase_and_Energy_for_Musical_Onset_Detection_in_the_Complex_Domain/links/5412b6110cf2bb7347dafd25/On-the-Use-of-Phase-and-Energy-for-Musical-Onset-Detection-in-the-Complex-Domain.pdf)
- [Estimation of Frequency, Amplitude, and Phase from the DFT of a Time Series](https://pdfs.semanticscholar.org/df2e/2b3ae9d784e19ea0840f8bb26ff622b17c22.pdf)
- [Positive Unlabeled Learning with Non-Negative Risk Estimator](http://papers.nips.cc/paper/6765-positive-unlabeled-learning-with-non-negative-risk-estimator.pdf)
