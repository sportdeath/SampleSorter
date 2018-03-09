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
    x, y = np.ogrid[:batch_size, :dataset.shape[1]]
    y = y - rotations[:,np.newaxis]
    batch = batch[x, y]

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

def octave_classifier_sum(batch_size, octave_length, reuse=False, training=False):
    # Octave classifier
    batch_positive_ph, decision_positive = octave_classifier(batch_size, octave_length, training=training, reuse=reuse)
    batch_negative_ph, decision_negative = octave_classifier(batch_size, octave_length, training=training, reuse=True)

    # Compute losses
    pi_p = 0.5
    beta = 0.01
    loss_positive = tf.reduce_mean(tf.sigmoid(-decision_positive))
    loss_unlabeled = tf.reduce_mean(tf.sigmoid(decision_negative))
    loss = pi_p * loss_positive + tf.maximum(-beta, loss_unlabeled - pi_p * (1 - loss_positive))

    # print(tf.get_collection(tf.GraphKeys.REGULARIZATION_LOSSES))
    # loss_regularizer = tf.reduce_mean(tf.get_collection(tf.GraphKeys.REGULARIZATION_LOSSES))
    # loss_regularizer_sum = tf.summary.scalar('loss_regularizer', loss_regularizer)

    # Make summaries
    loss_positive_sum = tf.summary.scalar('loss_positive', loss_positive)
    loss_unlabeled_sum = tf.summary.scalar('loss_unlabeled', loss_unlabeled)
    loss_sum = tf.summary.scalar('loss', loss)

    # Merge the summaries
    summary = tf.summary.merge((loss_positive_sum, loss_unlabeled_sum, loss_sum))

    # Optimize
    optimizer = None
    if training:
        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
        with tf.control_dependencies(update_ops):
            optimizer = tf.train.AdamOptimizer(0.0001).minimize(loss)

    return summary, optimizer, batch_positive_ph, batch_negative_ph

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
    num_train = int(num_octaves * 0.8)
    train_set = octaves[:num_train]
    validation_set = octaves[num_train:]

    train_summary, optimizer, batch_positive_ph, batch_negative_ph = octave_classifier_sum(
            batch_size, octave_length, training=True, reuse=False)
    validation_summary, _, batch_positive_val_ph, batch_negative_val_ph = octave_classifier_sum(
            batch_size, octave_length, training=False, reuse=True)

    with tf.Session() as session:
        session.run(tf.local_variables_initializer())
        session.run(tf.global_variables_initializer())

        train_writer = tf.summary.FileWriter("tmp/sort/11/train")
        validation_writer = tf.summary.FileWriter("tmp/sort/11/validation")
        train_writer.add_graph(session.graph)

        for i in range(100000000):
            # Make random positive and negative batches
            batch_positive = makePositiveBatch(train_set, batch_size)
            batch_negative = makeNegativeBatch(train_set, batch_size)
            batch_positive_val = makePositiveBatch(validation_set, batch_size)
            batch_negative_val = makeNegativeBatch(validation_set, batch_size)

            feed_dict = {}
            feed_dict[batch_positive_ph] = batch_positive
            feed_dict[batch_negative_ph] = batch_negative
            feed_dict[batch_positive_val_ph] = batch_positive_val
            feed_dict[batch_negative_val_ph] = batch_negative_val

            session.run(optimizer, feed_dict=feed_dict)

            if i % 100 == 0:
                train_summary_, validation_summary_ = session.run((train_summary, validation_summary), feed_dict=feed_dict)
                train_writer.add_summary(train_summary_, i)
                validation_writer.add_summary(validation_summary_, i)
                print(i)

if __name__ == "__main__":
    train()
