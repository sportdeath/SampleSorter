#! /usr/local/bin/python3

from sys import argv
import os

import numpy as np
import tensorflow as tf

from sample_reader import SampleReader

TENSORFLOW_GRAPH = "graph"
TUNING_EPSILON = 15
MOST_SEMITONES = 3

USER_LIBRARY = "A:\\User Library"
SAMPLE_LIBRARY = "A:\\User Library\\SampleLibrary\\"
FORCE_REPROCESS = False

def filter_samples(octaves, tempos, tunings, paths, sample_index):
    # Determine if tuning is ~semitone
    delta_tunings = 1200. * np.log2(tempos[sample_index]/tempos)
    delta_tunings = np.mod(delta_tunings + 600., 1200.) - 600.
    is_in_tune = (50. - np.abs(50. - np.mod(delta_tunings, 100.))) < TUNING_EPSILON
    is_in_tune = np.logical_and(is_in_tune, (np.abs(delta_tunings) < 100 * MOST_SEMITONES))

    octaves = octaves[is_in_tune]
    tempos = tempos[is_in_tune]
    tunings = tunings[is_in_tune]
    delta_tunings = delta_tunings[is_in_tune]
    paths = [path for (path, b) in zip(paths, is_in_tune) if b]

    # Roll the octave
    discrete_tunings = np.round(delta_tunings/100).astype(int)

    rows, cols = np.ogrid[:octaves.shape[0], :octaves.shape[1]]

    discrete_tunings = np.mod(discrete_tunings, octaves.shape[1])
    cols = cols - discrete_tunings[:,np.newaxis]

    octaves = octaves[rows, cols]

    tunings = tunings + delta_tunings

    return octaves, tempos, tunings, paths

def test_filter_samples():
    tempos = np.array([60, 67, 69, 59, 33.5, 90, 180])
    batch_size = tempos.shape[0]
    octaves = np.expand_dims(np.arange(12), axis=0)
    octaves = np.tile(octaves, (batch_size, 1))
    tunings = np.array(range(batch_size))
    paths = ['a'] * batch_size

    sample_index = 0

    octaves, tempos, tunings, paths = filter_samples(octaves, tempos, tunings, paths, sample_index)

def main():

    if len(argv) < 2:
        print("Input the sample name")
        return

    sample_path = argv[1]
    sample_path = os.path.abspath(sample_path)

    print("Loading octaves...")
    octaves, tempos, tunings, paths = SampleReader.readSamples(
            USER_LIBRARY,
            SAMPLE_LIBRARY,
            FORCE_REPROCESS,
            shuffle=False)
    print("Found " + str(octaves.shape[0]) + " octaves.")



    sample_index = paths.index(sample_path)
    sample_octave = octaves[sample_index]

    octaves, tempos, tunings, paths = filter_samples(octaves, tempos, tunings, paths, sample_index)

    saver = tf.train.import_meta_graph(TENSORFLOW_GRAPH + ".meta")

    with tf.Session() as session:
        saver.restore(session, TENSORFLOW_GRAPH)

        graph = tf.get_default_graph()
        octave_ph = graph.get_tensor_by_name("octave_classifier/octave:0")
        training_ph = graph.get_tensor_by_name("octave_classifier/training:0")
        decision = graph.get_tensor_by_name("octave_classifier/decision:0")

        decisions = []

        # for index, octave in enumerate(octaves):

        test_octaves = sample_octave + octaves
        test_octaves = test_octaves/np.linalg.norm(test_octaves, axis=1)[:, None]

        feed_dict = {}
        feed_dict[octave_ph] = test_octaves
        feed_dict[training_ph] = False
        decision_ = session.run(tf.sigmoid(decision), feed_dict=feed_dict)
        decision_ = list(decision_[:,0])

        zipped = zip(decision_, tunings, paths)
        sort = sorted(zipped, key=lambda x: x[0])

        for i in range(1,200):
            print(sort[-i])

if __name__ == "__main__":
    main()
