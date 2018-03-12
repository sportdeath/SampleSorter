# SampleSorter

This code analyzes a library of audio samples and determines each sample's tempo, tuning, and harmonic profile.
Using these values we determine pairs of samples that, when played at particular speeds, 
have the same tempo, have the same tuning and are harmonically coherent.
It is intended as a compositional tool in making sample-based music.

## Motivation

Musicians began making compositions out of other recordings as early as the 40's with the development of [musique concrete](https://en.wikipedia.org/wiki/Musique_concr%C3%A8te).
Sampling was primarily an experimental form of music for many decades, 
making rare appearances on pop records such as The Beatle's [Revolution 9](https://www.youtube.com/watch?v=HWmvbxGpra4).
Sampling broke into the mainstream the 80's with the advent of hip-hop and has become a staple of popular music ever since.

In the early days of hip-hop one of the dominant sampling techniques was to combine a [drum break and a loop](https://www.youtube.com/watch?v=q7Ej8Te_35g).
With more advanced sampling technology like the MPC-3000 and eventually the computer and digital workstations many more styles developed.
For example hip-hop pioneer
[J-Dilla](https://www.youtube.com/watch?v=hXeywtmWKzU)
developed a style that involved [chopping](https://en.wikipedia.org/wiki/Chopping_(sampling_technique)) up samples and rearranging the pieces.

The sampling technique this work serves to aid is sometimes known as "[plunderphonics](https://en.wikipedia.org/wiki/Plunderphonics)" or a "sound collage".
It was popularized by The Avalanches in their masterpiece [Since I Left You](https://www.youtube.com/watch?v=LhBacKKEyBU).
This technique involves layering many samples on top of each other.
DJs employ a related technique called [harmonic mixing](https://en.wikipedia.org/wiki/Harmonic_mixing) to transition smoothly between different songs.
Even though the samples are often untouched (no chopping or added effects), 
the combination of different moods and genres and the context they are placed in creates an entirely new sound.
There is something incredibly moving and surreal in being able to take recordings from different genres, different decades, and different continents and realizing them as different parts of the same song.

Despite its beauty, music made with this technique is rare, perhaps because of how difficult it is to create;
the process of finding groups of samples that can be played simultaneously is extremely time consuming.
This is partly because two fundamental musical properties, tempo and pitch, are inherently linked.
Lowering or raising the pitch of a song slows down or speeds up the tempo.
If two samples do in fact have the same tempo and tuning 
(tuning measured in cents as the deviation from [12-tone equal tempermant](https://en.wikipedia.org/wiki/Equal_temperament) with [A4 = 440hz](https://en.wikipedia.org/wiki/A440_(pitch_standard))),
they might not be harmonically coherent.
The samples could be in conflicting keys or modes, again greatly reducing the number valid sample pairs.
Once a musician has finally found samples that have consistent rhythm and harmony,
they then get to choose which samples should actually go together based on aesthetics and composition.

This code works to make this process easier by automatically detecting which samples could be repitched to have the same tempo and tuning. Each of these pairs is then given a rating of how "harmonically coherent" it is using a classifier.

It should be noted that a process called [time stretching](https://en.wikipedia.org/wiki/Audio_time_stretching_and_pitch_scaling) is able to change the speed of audio without affecting the pitch.
However it necessitates windowing the audio which introduces digital artifacts particularly at transients.
Modern algorithms cover up most of these artifacts, but there is nevertheless something uncanny about time-stretched audio even for small changes in speed.

## Problem Setup

I have a collection of audio samples in Ableton's (\*.alc file format)[https://help.ableton.com/hc/en-us/articles/209769625-Live-specific-file-types-adg-als-alp-] (a convenient way to make references to sections of audio).
These are small (~3-30s long) snippets of music from a variety of genres; 
soul, jazz, funk, rock, hip-hop, disco, electronic, experimental, etc.
My own collection is not public due to copyright, 
but there is a large online database of the samples used in popular music at [WhoSampled.com](https://www.whosampled.com/most-sampled-tracks/)

I want to be able to layer these samples to create new music.
To layer sounds I can do any of the following:

- Adjust the speed at which the samples are played. This affects both the pitch and tempo.
- Choose when each sample begins to play.
- Adjust the volume of each sample.
- Add filters and effects to each sample.

Any pair of samples I choose to layer must abide by the following:

- The samples need to have the same tempo.
- The samples need to be in tune.
- The combination must sound "musical".

For the purpose of this project, I only care about adjusting the relative speed and offset between samples.
Mixing, filtering and adding effects is left to the musician.

## Solution Description

### Tempo

Autocorrelation

### Octave-Mapped Spectrogram

#### Tuning

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
