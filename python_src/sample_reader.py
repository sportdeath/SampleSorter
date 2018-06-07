import os
import tensorflow as tf
import numpy as np

from sample import Sample
from octave_classifier import OctaveClassifier

class SampleReader:
    @staticmethod
    def read_samples(user_library, sample_library, force_reprocess, shuffle=True):
        user_library = os.path.expanduser(user_library)
        sample_library = os.path.expanduser(sample_library)

        octaves = []
        tempos = []
        tunings = []
        paths = []
        audio_file_paths = []
        for dirpath, dirnames, filenames in os.walk(sample_library):
            for filename in filenames:
                name, ext = os.path.splitext(filename)
                if ext == ".alc":
                    full_path = os.path.abspath(os.path.join(dirpath, filename))
                    filepath = os.path.join(dirpath, filename)
                    octave, tempo, tuning, audio_file_path = SampleReader.read_sample(filepath, user_library, force_reprocess)
                    if np.count_nonzero(octave) == 0 or np.any(np.isnan(octave)):
                        print("Cannot read ", filepath, "\n")
                        continue
                    paths.append(full_path)
                    octaves.append(octave)
                    tempos.append(tempo)
                    tunings.append(tuning)
                    audio_file_paths.append(audio_file_path)

        if shuffle:
            zipped = list(zip(octaves, tempos, tunings, paths, audio_file_paths))
            np.random.shuffle(zipped)
            octaves, tempos, tunings, paths, audio_file_paths = zip(*zipped)

        octaves = np.array(octaves)
        tempos = np.array(tempos)
        tunings = np.array(tunings)

        # Normalize
        norm = np.linalg.norm(octaves, axis=1)
        octaves = octaves/norm[:,None]

        return octaves, tempos, tunings, paths, audio_file_paths

    @staticmethod
    def read_sample(filepath, user_library, force_reprocess):
        s = Sample(filepath, user_library, force_reprocess)
        s.process()
        octave = s.getOctave()
        tempo = s.getTempo()
        tuning = s.getTuning()
        audio_path = s.getAudioPath()
        # if force_reprocess:
        s.writeToFile()
        s.delete()
        return octave, tempo, tuning, audio_path

    @staticmethod
    def make_positive_batch(dataset, batch_size):
        random_sel = np.random.randint(dataset.shape[0], size=batch_size)
        batch = dataset[random_sel]

        # Rotate each by a random amount
        rotations = np.random.randint(12, size=batch_size)
        x, y = np.ogrid[:batch_size, :dataset.shape[1]]
        y = y - rotations[:,np.newaxis]
        batch = batch[x, y]

        return batch

    @staticmethod
    def make_unlabeled_batch(dataset, batch_size):
        batch0 = SampleReader.make_positive_batch(dataset, batch_size)
        batch1 = SampleReader.make_positive_batch(dataset, batch_size)

        # Add them together
        batch = batch0 + batch1

        # normalize
        norm = np.linalg.norm(batch, axis=1)
        batch = batch/norm[:,None]

        return batch

    @staticmethod
    def make_batch(dataset, batch_size):
        positive_batch = SampleReader.make_positive_batch(dataset, batch_size)
        unlabeled_batch = SampleReader.make_unlabeled_batch(dataset, batch_size)

        batch = np.concatenate((positive_batch, unlabeled_batch))

        return batch
