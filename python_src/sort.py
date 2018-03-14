import numpy as np
import tensorflow as tf

class OctaveClassifier:
    def __init__(
            self,
            activation=tf.nn.relu,
            octave_length=12,
            units_per_layer=100,
            localization_layers=3,
            classification_layers=5,
            dropout=0.8,
            expected_positive=0.5,
            non_negative_tolerance=0.,
            learning_rate=0.0005):

        self.activation = activation
        self.octave_length = octave_length
        self.localization_layers = [units_per_layer] * localization_layers + [self.octave_length]
        self.classification_layers = [units_per_layer] * classification_layers + [1]
        self.dropout = dropout
        self.expected_positive = expected_positive
        self.non_negative_tolerance = non_negative_tolerance
        self.learning_rate = learning_rate

    def construct_with_loss(self, batch_size, reuse=False):
        # Make training placeholder
        training = tf.placeholder(tf.bool)

        # Construct classifier
        batch_ph, decision, summaries = self.construct(batch_size, training=training, reuse=reuse)

        summaries = tf.summary.merge(summaries)

        # Compute losses
        loss, loss_summaries = self.loss(decision)

        # Optimize
        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS) # For batch norm
        with tf.control_dependencies(update_ops):
            optimizer = tf.train.AdamOptimizer(self.learning_rate).minimize(loss)

        return optimizer, loss_summaries, batch_ph, training, summaries

    def loss(self, decision, name="octave_classifier_loss"):
        """
        Compute the non-negative loss as described in
        https://arxiv.org/abs/1703.00593
        """
        with tf.variable_scope(name):
            batch_size = decision.get_shape().as_list()[0]
            positive_size = int(batch_size/2)

            decision_positive = decision[:positive_size]
            decision_unlabeled = decision[positive_size:]

            # Regress the positive and unlabeled decisions
            loss_positive = tf.reduce_mean(tf.sigmoid(-decision_positive))
            loss_unlabeled = tf.reduce_mean(tf.sigmoid(decision_unlabeled))

            # Compute the negative loss
            loss_negative = loss_unlabeled - self.expected_positive * (1 - loss_positive)

            # Take the maximum so the negative loss is non negative
            loss_non_negative = tf.maximum(-self.non_negative_tolerance, loss_negative)

            # Sum the losses
            loss = (self.expected_positive * loss_positive) + loss_non_negative

            # Make summaries
            loss_positive_summary = tf.summary.scalar('loss_positive', loss_positive)
            loss_unlabeled_summary = tf.summary.scalar('loss_unlabeled', loss_unlabeled)
            loss_negative_summary = tf.summary.scalar('loss_negative', loss_negative)
            loss_non_negative_summary = tf.summary.scalar('loss_non_negative', loss_non_negative)
            loss_summary = tf.summary.scalar('loss', loss)

            # Merge the summaries
            loss_summaries = tf.summary.merge((
                loss_positive_summary, 
                loss_unlabeled_summary, 
                loss_negative_summary,
                loss_non_negative_summary,
                loss_summary))

        return loss, loss_summaries

    def construct(self, batch_size, training, name="octave_classifier", reuse=False):
        with tf.variable_scope(name, reuse=reuse):
            octave = tf.placeholder(dtype=tf.float32, shape=[2 * batch_size, self.octave_length])

            # Determine how the octave should be rotated
            choices, summaries0 = self._dense_net(octave, self.localization_layers, training, "localization", reuse=reuse)
            
            # Rotate the octave
            transformed = self._octave_rotate_disc(octave, choices)

            # Classify the octave
            decision, summaries1 = self._dense_net(transformed, self.classification_layers, training, "classification", reuse=reuse)

            summaries = summaries0 + summaries1

        return octave, decision, summaries

    def _dense_net(self, input_, layer_units, training, name="dense_net", reuse=False):
        """
        Make a dense neural net where each layer is an entry in
        layer_units. All but the last layer includes a nonlinearity.
        """
        hidden = input_

        with tf.variable_scope(name, reuse=reuse):

            summaries = []

            for i, num_units in enumerate(layer_units):
                # Make the last layer linear
                activation = None
                if i < len(layer_units) - 1:
                    activation = self.activation

                # Dense connection
                hidden = tf.layers.dense(
                        inputs=hidden,
                        units=num_units,
                        activation=activation,
                        kernel_initializer=tf.contrib.layers.xavier_initializer(),
                        name="dense_" + str(i),
                        reuse=reuse)

                mean, variance = tf.nn.moments(hidden, [0])
                mean_summary = tf.summary.histogram('mean_' + str(i), mean)
                variance_summary = tf.summary.histogram('variance_' + str(i), variance)
                summaries.append(mean_summary)
                summaries.append(variance_summary)

                if i < len(layer_units) - 1:
                    # Batch renorm
                    # https://arxiv.org/pdf/1702.03275.pdf
                    hidden = tf.layers.batch_normalization(
                            hidden, 
                            training=training, 
                            name="bn_" + str(i), 
                            renorm=True,
                            fused=True,
                            reuse=reuse)

                    dropout = tf.where(training, self.dropout, 1)
                    hidden = tf.nn.dropout(hidden, dropout)

        return hidden, summaries

    def _octave_rotate_disc(self, octave, choices, name="octave_rotate_disc"):
        with tf.variable_scope(name):
            # Extract the shape
            octave_shape = octave.get_shape().as_list()
            batch_size = octave_shape[0]
            octave_length = octave_shape[1]

            # Make the indices that will rotate each octave
            # by each possible integer
            indices = []
            for i in range(octave_length):
                indices.append(np.roll(np.arange(octave_length), -i))
            indices = np.array(indices)
            # [ 0  1 ... 11 12]
            # [ 1  2 ... 11  0]
            #        ...
            # [12  0 ... 10 11]
            indices = np.expand_dims(indices, axis=0)
            indices = np.tile(indices, (batch_size, 1, 1))

            batch = np.arange(batch_size)
            batch = np.expand_dims(np.expand_dims(batch, axis=-1), axis=-1)
            batch = np.tile(batch, (1, octave_length, octave_length))

            # Stack them
            indices = np.stack((batch, indices), axis=-1)

            # Convert to tf
            indices = tf.constant(indices, tf.int32)

            # Rotate each octave
            octave_rotations = tf.gather_nd(octave, indices)

            # Make the choices
            choices = tf.nn.softmax(choices)
            choices = tf.expand_dims(choices, axis=-1)
            choices = tf.tile(choices, (1, 1, octave_length))
            octave_rotated = octave_rotations * choices

            # Sum the choices
            octave_rotated = tf.reduce_sum(octave_rotated, axis=1)

            return octave_rotated

    def _test_octave_rotate_disc(self):
        octave = tf.constant(
                [[1, 2, 3, 4, 5],
                 [5, 4, 3, 2, 1]], tf.float32)

        choices = tf.constant(
                [[10., 1., 2., 5., 1.],
                 [1.0, 7., 1., 20., 1.]], tf.float32)

        _octave_rotate_disc(octave, choices)

if __name__=="__main__":
    _test_octave_rotate_disc()
