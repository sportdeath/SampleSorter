import numpy as np
import tensorflow as tf

def octave_classifier(batch_size, octave_length, name="octave_classifier", reuse=False):
    with tf.variable_scope(name):
        octave = tf.placeholder(dtype=tf.float32, shape=[batch_size, octave_length])

        localization_layers = [100, 100, 100, 2]
        xy = _dense_net(octave, localization_layers, "localization", reuse=reuse)

        with tf.variable_scope("octave_transform"):
            # Compute the angle and scale from the localization
            angle = tf.atan2(xy[:,0], xy[:,1])
            scale = tf.norm(xy, axis=1)
            # Make it [batch_size, 1]
            angle = tf.expand_dims(angle, axis=1)
            scale = tf.expand_dims(scale, axis=1)

            # Rotate the octave by the angle
            transformed = scale * _octave_rotate(octave, angle)

        classification_layers = [100, 100, 100, 100, 100, 1]
        decision = _dense_net(transformed, classification_layers, "classification", reuse=reuse)
        decision = tf.sigmoid(decision)

    return octave, decision, scale

def _dense_net(input_, layer_units, name="dense_net", reuse=False):
    """
    Make a dense neural net where each layer is an entry in
    layer_units. All but the last layer includes a nonlinearity.
    """
    hidden = input_

    with tf.variable_scope(name):

        for i, num_units in enumerate(layer_units):
            # Make the last layer linear
            activation = None
            if i < len(layer_units) - 1:
                activation = tf.nn.relu

            # Dense connection
            hidden = tf.layers.dense(
                    inputs=hidden,
                    units=num_units,
                    activation=activation,
                    kernel_initializer=tf.contrib.layers.xavier_initializer(),
                    name="dense_" + str(i),
                    reuse=reuse)

    return hidden

def _octave_rotate(octave, angle, name="octave_rotate"):
    """
    Bilinearly rotate the octave by the angle.
    """

    with tf.variable_scope(name):
        # Extract the shape
        octave_shape = octave.get_shape().as_list()
        batch_size = octave_shape[0]
        octave_length = octave_shape[1]

        # Compute the displacement from the angle
        displacement = (angle * octave_length)/(2 * np.pi)

        # Establish basic indices
        indices = tf.range(octave_length, dtype=tf.float32) + displacement
        indices = tf.mod(indices, octave_length)

        # Add the diplacement to the indices
        indices_f = tf.floor(indices)
        indices_c = tf.mod(indices_f + 1, octave_length)
        # Compute bilinear interpolation weights
        weight_c = indices - indices_f
        weight_f = 1 - weight_c
        # Convert to int
        indices_f = tf.cast(indices_f, tf.int32)
        indices_c = tf.cast(indices_c, tf.int32)

        # Add batch to the indices
        batch = tf.expand_dims(tf.range(batch_size, dtype=tf.int32), axis=-1)
        batch = tf.tile(batch, [1, octave_length])
        indices_f = tf.stack((batch, indices_f), axis=-1)
        indices_c = tf.stack((batch, indices_c), axis=-1)

        # Displace and interpolate
        octave_out = weight_f * tf.gather_nd(octave, indices_f) + \
                     weight_c * tf.gather_nd(octave, indices_c)

    return octave_out

def _test_octave_rotate(batch_size, octave_length):
    # Each octave in the batch is [0..11]
    octave = tf.expand_dims(tf.range(octave_length, dtype=tf.float32), 0)
    octave = tf.tile(octave, [batch_size, 1])
    print("Octave:\n", tf.Session().run(octave))

    # The angles are [[0...2Pi]]
    angle = tf.range(batch_size, dtype=tf.float32)
    angle = angle * (2 * np.pi)/batch_size
    angle = tf.expand_dims(angle, axis=-1)
    print("Angle:\n", tf.Session().run(angle))

    out = _octave_rotate(octave, angle)
    print("Rotation:\n", tf.Session().run(out))
