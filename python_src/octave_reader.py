import os
import tensorflow as tf
import numpy as np

from sample import Sample
from octave_classifier import OctaveClassifier

class OctaveReader:

    @staticmethod
    def getOctave(filepath, user_library, force_reprocess):
        s = Sample(filepath, user_library, force_reprocess)
        s.process()
        octave = s.getOctave()
        s.writeToFile()
        s.delete()
        return octave

    @staticmethod
    def getOctaves(user_library, sample_library, force_reprocess):
        user_library = os.path.expanduser(user_library)
        sample_library = os.path.expanduser(sample_library)

        octaves = []
        for dirpath, dirnames, filenames in os.walk(sample_library):
            for filename in filenames:
                name, ext = os.path.splitext(filename)
                if ext == ".alc":
                    filepath = os.path.join(dirpath, filename)
                    octave = OctaveReader.getOctave(filepath, user_library, force_reprocess)
                    if np.count_nonzero(octave) == 0 or np.any(np.isnan(octave)):
                        print("Cannot read ", filepath, "\n")
                        continue
                    octaves.append(octave)

        # Convert to a numpy array and shuffle
        octaves = np.array(octaves)
        np.random.shuffle(octaves)

        # Normalize
        norm = np.linalg.norm(octaves, axis=1)
        octaves = octaves/norm[:,None]

        return octaves

    @staticmethod
    def makePositiveBatch(dataset, batch_size):
        random_sel = np.random.randint(dataset.shape[0], size=batch_size)
        batch = dataset[random_sel]

        # Rotate each by a random amount
        rotations = np.random.randint(12, size=batch_size)
        x, y = np.ogrid[:batch_size, :dataset.shape[1]]
        y = y - rotations[:,np.newaxis]
        batch = batch[x, y]

        return batch

    @staticmethod
    def makeUnlabeledBatch(dataset, batch_size):
        batch0 = OctaveReader.makePositiveBatch(dataset, batch_size)
        batch1 = OctaveReader.makePositiveBatch(dataset, batch_size)

        # Add them together
        batch = batch0 + batch1

        # normalize
        norm = np.linalg.norm(batch, axis=1)
        batch = batch/norm[:,None]

        return batch

    @staticmethod
    def makeBatch(dataset, batch_size):
        positive_batch = OctaveReader.makePositiveBatch(dataset, batch_size)
        unlabeled_batch = OctaveReader.makeUnlabeledBatch(dataset, batch_size)

        batch = np.concatenate((positive_batch, unlabeled_batch))

        return batch
