import os
import tensorflow as tf
import numpy as np

from sample import Sample
from sort import OctaveClassifier

def getOctave(filepath, user_library, force_reprocess):
    s = Sample(filepath, user_library, force_reprocess)
    s.process()
    octave = s.getOctave()
    s.writeToFile()
    s.delete()
    return octave

def getOctaves(user_library, sample_library, force_reprocess):
    user_library = os.path.expanduser(user_library)
    sample_library = os.path.expanduser(sample_library)

    octaves = []
    for dirpath, dirnames, filenames in os.walk(sample_library):
        for filename in filenames:
            name, ext = os.path.splitext(filename)
            if ext == ".alc":
                filepath = os.path.join(dirpath, filename)
                octave = getOctave(filepath, user_library, force_reprocess)
                if np.count_nonzero(octave) == 0:
                    print(filepath)
                    print(octave)
                octaves.append(octave)

    # Convert to a numpy array and shuffle
    octaves = np.array(octaves)
    np.random.shuffle(octaves)

    # mean = np.expand_dims(octaves.mean(axis=1), axis=1)
    # std = np.expand_dims(octaves.std(axis=1), axis=1)
    # octaves = (octaves - mean)/std
    # octaves = 
    norm = np.linalg.norm(octaves, axis=1)
    octaves = octaves/norm[:,None]

    return octaves

def makePositiveBatch(dataset, batch_size):
    random_sel = np.random.randint(dataset.shape[0], size=batch_size)
    batch = dataset[random_sel]

    # Rotate each by a random amount
    rotations = np.random.randint(12, size=batch_size)
    x, y = np.ogrid[:batch_size, :dataset.shape[1]]
    y = y - rotations[:,np.newaxis]
    batch = batch[x, y]

    return batch

def makeUnlabeledBatch(dataset, batch_size):
    batch0 = makePositiveBatch(dataset, batch_size)
    batch1 = makePositiveBatch(dataset, batch_size)

    # Add them together
    batch = batch0 + batch1

    # normalize
    norm = np.linalg.norm(batch, axis=1)
    batch = batch/norm[:,None]

    return batch

def train():
    batch_size = 5
    force_reprocess = False
    user_library = "~/Documents/UserLibrary/"
    sample_library = "~/Documents/UserLibrary/SampleLibrary/"

    # Get all the octaves
    octaves = getOctaves(user_library, sample_library, force_reprocess)

    # Get the length (12)
    octave_length = octaves.shape[1]

    # Get training and validation sets
    num_octaves = octaves.shape[0]
    num_train = int(num_octaves * 0.8)
    train_set = octaves[:num_train]
    print(train_set.shape)
    validation_set = octaves[num_train:]
    print(validation_set.shape)

    oc = OctaveClassifier()
    optimizer, summary, batch_positive_ph, batch_unlabeled_ph, training = oc.train(batch_size)

    with tf.Session() as session:
        session.run(tf.local_variables_initializer())
        session.run(tf.global_variables_initializer())

        train_writer = tf.summary.FileWriter("tmp/sort/32/train")
        validation_writer = tf.summary.FileWriter("tmp/sort/32/validation")
        train_writer.add_graph(session.graph)

        saver = tf.train.Saver()

        for i in range(100000000):
            # Make random positive and unlabeled batches
            batch_positive = makePositiveBatch(train_set, batch_size)
            batch_unlabeled = makeUnlabeledBatch(train_set, batch_size)

            feed_dict = {}
            feed_dict[batch_positive_ph] = batch_positive
            feed_dict[batch_unlabeled_ph] = batch_unlabeled
            feed_dict[training] = True
            _, train_summary = session.run((optimizer, summary), feed_dict=feed_dict)

            if i % 100 == 0:
                feed_dict={}
                batch_positive_val = makePositiveBatch(validation_set, batch_size)
                batch_unlabeled_val = makeUnlabeledBatch(validation_set, batch_size)
                feed_dict[batch_positive_ph] = batch_positive_val
                feed_dict[batch_unlabeled_ph] = batch_unlabeled_val
                feed_dict[training] = False

                validation_summary = session.run(summary, feed_dict=feed_dict)

                train_writer.add_summary(train_summary, i)
                validation_writer.add_summary(validation_summary, i)
                print(i)

            if i % 100000 == 0:
                print("Writing...")
                saver.save(session, "tmp/sort/32/model.ckpt")

if __name__ == "__main__":
    train()
