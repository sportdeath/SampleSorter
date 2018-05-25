import itertools
import tensorflow as tf

from octave_classifier import OctaveClassifier
from octave_reader import OctaveReader

USER_LIBRARY = "~/FatDisk/User Library/"
SAMPLE_LIBRARY = "~/Samples/"
FORCE_REPROCESS = False
TRAIN_PROPORTION = 0.8
LOG_DIR = "tmp/octave_classifier/hyper_param_search2/"

HYPERPARAMS = [
        [10], # batch_size
        [0.0005], # learning_rate
        [0.5], # dropout
        [20, 50], # units_per_layer
        [2], # localization_layers
        [2, 3], # classification_layers
        ]

def make_hyperparam_string(hyperparams):
    s = ""
    s += "bn=" + str(hyperparams[0])
    s += "lr=" + str(hyperparams[1])
    s += "drop=" + str(hyperparams[2])
    s += "upl=" + str(hyperparams[3])
    s += "loc=" + str(hyperparams[4])
    s += "clas=" + str(hyperparams[5])
    return s

def main():
    # Get all the octaves
    print("Loading octaves...")
    octaves = OctaveReader.getOctaves(USER_LIBRARY, SAMPLE_LIBRARY, FORCE_REPROCESS)
    print("Done.")

    # Get training and validation sets
    num_train = int(octaves.shape[0] * TRAIN_PROPORTION)
    train_set = octaves[:num_train]
    validation_set = octaves[num_train:]

    # For each set of hyperparameters
    for hyperparams in itertools.product(*HYPERPARAMS):
        print("Using params:")
        print(make_hyperparam_string(hyperparams))

        batch_size = hyperparams[0]

        # Reset the graph
        tf.reset_default_graph()

        # Generate the graph
        oc = OctaveClassifier(
                learning_rate=hyperparams[1],
                dropout=hyperparams[2],
                units_per_layer=hyperparams[3],
                localization_layers=hyperparams[4],
                classification_layers=hyperparams[5]
                )
        optimizer, summary, batch_ph, training = oc.construct_with_loss(batch_size)

        log_dir = LOG_DIR + make_hyperparam_string(hyperparams) + "/"

        with tf.Session() as session:
            session.run(tf.local_variables_initializer())
            session.run(tf.global_variables_initializer())

            train_writer = tf.summary.FileWriter(log_dir + "train")
            validation_writer = tf.summary.FileWriter(log_dir + "validation")

            for i in range(3000000):
                # Make random positive and unlabeled batches
                batch = OctaveReader.makeBatch(train_set, batch_size)

                feed_dict = {}
                feed_dict[batch_ph] = OctaveReader.makeBatch(train_set, batch_size)
                feed_dict[training] = True

                session.run(optimizer, feed_dict=feed_dict)

                if i % 2000 == 0:
                    train_summary = session.run(summary, feed_dict=feed_dict)

                    feed_dict={}
                    feed_dict[batch_ph] = OctaveReader.makeBatch(validation_set, batch_size)
                    feed_dict[training] = False

                    validation_summary = session.run(summary, feed_dict=feed_dict)

                    train_writer.add_summary(train_summary, i)
                    validation_writer.add_summary(validation_summary, i)
                    print(i)

if __name__ == "__main__":
    main()
