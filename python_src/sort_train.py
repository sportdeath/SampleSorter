import os
import tensorflow as tf
import numpy as np

from sample import Sample
from sort import octave_classifier

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

    return octaves

def makePositiveBatch(dataset, batch_size):
    random_sel = np.random.randint(dataset.shape[0], size=batch_size)
    batch = dataset[random_sel]

    # Normalize each one
    norm = np.linalg.norm(batch, axis=1)
    batch = batch/norm[:,None]

    # Rotate each by a random amount
    rotations = np.random.randint(12, size=batch_size)
    batch = np.roll(batch, rotations, axis=1)

    return batch

def makeNegativeBatch(dataset, batch_size):
    batch0 = makePositiveBatch(dataset, batch_size)
    batch1 = makePositiveBatch(dataset, batch_size)

    # Add them together
    batch = batch0 + batch1

    # normalize
    norm = np.linalg.norm(batch, axis=1)
    batch = batch/norm[:,None]

    return batch

def train():
    batch_size = 10
    force_reprocess = False
    user_library = "~/Documents/UserLibrary/"
    sample_library = "~/Documents/UserLibrary/SampleLibrary/"

    # Get all the octaves
    octaves = getOctaves(user_library, sample_library, force_reprocess)

    # Get the length (12)
    octave_length = octaves.shape[1]

    # Get training and validation sets
    num_octaves = octaves.shape[0]
    num_train = int(num_octaves/2)
    train_set = octaves[:num_train]
    validation_set = octaves[num_train:]

    # Octave classifier
    batch_positive_ph, decision_positive, scale_positive = octave_classifier(batch_size, octave_length)
    batch_negative_ph, decision_negative, scale_negative = octave_classifier(batch_size, octave_length, reuse=True)

    # Compute losses
    loss_scale = tf.reduce_mean(tf.square(scale_positive - 1) + tf.square(scale_negative - 1))
    tf.summary.scalar('loss_scale', loss_scale)
    loss_positive = tf.reduce_mean(-tf.log(decision_positive))
    tf.summary.scalar('loss_positive', loss_positive)
    loss_negative = tf.reduce_mean(-tf.log(1 - decision_negative))
    tf.summary.scalar('loss_negative', loss_negative)
    loss = loss_positive + loss_negative + loss_scale
    tf.summary.scalar('loss', loss)

    # Optimize
    optimizer = tf.train.AdamOptimizer(0.0001).minimize(loss)

    merged_summary = tf.summary.merge_all()

    with tf.Session() as session:
        session.run(tf.local_variables_initializer())
        session.run(tf.global_variables_initializer())

        writer = tf.summary.FileWriter("tmp/sort/14")
        writer.add_graph(session.graph)

        for i in range(100000000):
            # Make random positive and negative batches
            batch_positive = makePositiveBatch(train_set, batch_size)
            batch_negative = makeNegativeBatch(train_set, batch_size)

            feed_dict = {}
            feed_dict[batch_positive_ph] = batch_positive
            feed_dict[batch_negative_ph] = batch_negative

            _, summary = session.run((optimizer, merged_summary), feed_dict=feed_dict)

            if i % 100 == 0:
                writer.add_summary(summary, i)
                print(i)

if __name__ == "__main__":
    train()
