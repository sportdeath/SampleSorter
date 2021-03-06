import numpy as np
import tensorflow as tf

class OctaveClassifier:
    def __init__(
            self,
            activation=tf.nn.relu,
            octave_length=12,
            units_per_layer=50,
            localization_layers=2,
            classification_layers=2,
            dropout=0.5,
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

    def construct_with_loss(self, reuse=False):

        # Construct classifier
        batch_ph, training_ph, decision = self.construct(reuse=reuse)

        # Compute losses
        loss, loss_summaries = self.loss(decision)

        # Optimize
        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS) # For batch norm
        with tf.control_dependencies(update_ops):
            optimizer = tf.train.AdamOptimizer(self.learning_rate).minimize(loss)

        return optimizer, loss_summaries, batch_ph, training_ph

    def loss(self, decision, name="octave_classifier_loss"):
        """
        Compute the non-negative loss as described in
        https://arxiv.org/abs/1703.00593
        """
        with tf.variable_scope(name):
            positive_size = tf.shape(decision)[0]/2
            positive_size = tf.cast(positive_size, tf.int32)

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

    def construct(self, name="octave_classifier", reuse=False):
        with tf.variable_scope(name, reuse=reuse):
            octave = tf.placeholder(dtype=tf.float32, name="octave")
            octave = tf.reshape(octave, [-1, self.octave_length])
            training = tf.placeholder(tf.bool, name="training")

            # Determine how the octave should be rotated
            choices = self._dense_net(octave, self.localization_layers, training, "localization", reuse=reuse)
            
            # Rotate the octave
            transformed = self._octave_rotate_disc(octave, choices)

            # Classify the octave
            decision = self._dense_net(transformed, self.classification_layers, training, "classification", reuse=reuse)
            decision = tf.identity(decision, name="decision")

        return octave, training, decision

    def _dense_net(self, input_, layer_units, training, name="dense_net", reuse=False):
        """
        Make a dense neural net where each layer is an entry in
        layer_units. All but the last layer includes a nonlinearity.
        """
        hidden = input_

        with tf.variable_scope(name, reuse=reuse):

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

                    # Dropout only if training
                    dropout = tf.where(training, self.dropout, 1)
                    hidden = tf.nn.dropout(hidden, dropout)

        return hidden

    def _octave_rotate_disc(self, octave, choices, name="octave_rotate_disc"):
        with tf.variable_scope(name):
            # Make the indices that will rotate each octave
            # by each possible integer
            indices = []
            for i in range(self.octave_length):
                indices.append(np.roll(np.arange(self.octave_length), -i))
            indices = np.array(indices)
            # [ 0  1 ... 11 12]
            # [ 1  2 ... 11  0]
            #        ...
            # [12  0 ... 10 11]
            indices = np.expand_dims(indices, axis=0)
            indices = tf.constant(indices, tf.int32)
            indices = tf.tile(indices, (tf.shape(octave)[0], 1, 1))

            batch = tf.range(tf.shape(octave)[0])
            batch = tf.expand_dims(tf.expand_dims(batch, axis=-1), axis=-1)
            batch = tf.tile(batch, (1, self.octave_length, self.octave_length))

            # Stack them
            indices = tf.stack((batch, indices), axis=-1)

            # Rotate each octave
            octave_rotations = tf.gather_nd(octave, indices)

            # Make the choices
            choices = tf.nn.softmax(choices)
            choices = tf.expand_dims(choices, axis=-1)
            choices = tf.tile(choices, (1, 1, self.octave_length))
            octave_rotated = octave_rotations * choices

            # Sum the choices
            octave_rotated = tf.reduce_sum(octave_rotated, axis=1)

            return octave_rotated

    def _test_octave_rotate_disc(self):
        octave = tf.constant(
                [[1, 2, 3, 4, 5],
                 [5, 4, 3, 2, 1]], tf.float32)
        octave = tf.reshape(octave, (-1, 5))

        choices = tf.constant(
                [[10., 1., 2., 5., 1.],
                 [1.0, 7., 1., 20., 1.]], tf.float32)
        choices = tf.reshape(choices, (-1, 5))

        self._octave_rotate_disc(octave, choices)

if __name__=="__main__":
    OctaveClassifier(octave_length=5)._test_octave_rotate_disc()
