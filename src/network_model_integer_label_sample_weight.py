# network_model_integer_label_sample_weight.py
#
# Name: Laura Tung
#
# Usage: python network_model_integer_label_sample_weight.py
#


#import pdb; pdb.set_trace() # Uncomment to debug code using pdb (like gdb)

import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle

from tensorflow import keras
from keras.models import Sequential, Model
from keras.layers import Input, Activation, BatchNormalization, Dropout, Dense, Bidirectional, LSTM, TimeDistributed, Concatenate
from keras.layers import Masking


def build_model(left_kmers_input_shape, right_kmers_input_shape, read_features_input_shape, consensus_features_input_shape, rest_features_input_shape):
    
    left_kmers_input = Input(shape=left_kmers_input_shape, dtype='float', name='left_kmers_input')
    right_kmers_input = Input(shape=right_kmers_input_shape, dtype='float', name='right_kmers_input')
    read_features_input = Input(shape=read_features_input_shape, dtype='float', name='read_features_input')
    consensus_features_input = Input(shape=consensus_features_input_shape, dtype='float', name='consensus_features_input')
    rest_features_input = Input(shape=rest_features_input_shape, dtype='float', name='rest_features_input')
    
    # Masking
    masking_left_kmers = Masking(mask_value=-1., input_shape=left_kmers_input_shape)(left_kmers_input)
    masking_right_kmers = Masking(mask_value=-1., input_shape=right_kmers_input_shape)(right_kmers_input)
    masking_read_features = Masking(mask_value=-1., input_shape=read_features_input_shape)(read_features_input)
    masking_consensus_features = Masking(mask_value=-1., input_shape=consensus_features_input_shape)(consensus_features_input)
    masking_rest_features = Masking(mask_value=-1., input_shape=rest_features_input_shape)(rest_features_input)
    
    # TimeDistributed Batch Normalization
    batch_norm_left_kmers = TimeDistributed(BatchNormalization())(masking_left_kmers)
    batch_norm_right_kmers = TimeDistributed(BatchNormalization())(masking_right_kmers)
    
    # TimeDistributed LSTM
    hier_lstm_left_kmers = TimeDistributed(LSTM(16, name='hier_lstm_left_kmers'))(batch_norm_left_kmers)
    hier_lstm_right_kmers = TimeDistributed(LSTM(16, name='hier_lstm_right_kmers'))(batch_norm_right_kmers)
    
    # Concatenate read_features, left_kmers, and right_kmers
    concat1 = Concatenate(axis=-1)([masking_read_features, hier_lstm_left_kmers, hier_lstm_right_kmers])
    
    # Batch Normalization
    batch_norm_read_branch = BatchNormalization()(concat1)
    batch_norm_consensus_branch = BatchNormalization()(masking_consensus_features)
    
    # Bi-LSTM
    bi_lstm_read_branch = Bidirectional(LSTM(64, return_sequences=True, dropout=0.1, name='bi_lstm_read_branch'))(batch_norm_read_branch)
    bi_lstm_consensus_branch = Bidirectional(LSTM(64, return_sequences=True, dropout=0.1, name='bi_lstm_consensus_branch'))(batch_norm_consensus_branch)
    
    # Concatenate read_branch, consensus_branch, and rest_features
    concat2 = Concatenate(axis=-1)([bi_lstm_read_branch, bi_lstm_consensus_branch, masking_rest_features])
    
    # TimeDistributed Dense 1
    dense1 = TimeDistributed(Dense(128, activation='relu', name='dense1'))(concat2)
    
    # TimeDistributed Dropout 1
    dropout1 = TimeDistributed(Dropout(.2))(dense1)
    
    # TimeDistributed Dense 2
    dense2 = TimeDistributed(Dense(32, activation='relu', name='dense2'))(dropout1)    

    # TimeDistributed Dropout 2
    dropout2 = TimeDistributed(Dropout(.2))(dense2)
    
    # TimeDistributed Dense softmax: final output
    output_dense = TimeDistributed(Dense(10, activation='softmax', name='output_dense'))(dropout2)
    
    model = Model(inputs=[left_kmers_input, right_kmers_input, read_features_input, consensus_features_input, rest_features_input], outputs=output_dense)
    
    print(model.summary())
    
    model.compile(optimizer='adam',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'],
                  sample_weight_mode="temporal")
    
    return model


def train_model(model, train_left_kmers_input, train_right_kmers_input, train_read_features_input, train_consensus_features_input, train_rest_features_input, train_true_labels, epochs, batch_size, validation_split, verbose, sample_weights):
    
    history = model.fit([train_left_kmers_input, train_right_kmers_input, train_read_features_input, train_consensus_features_input, train_rest_features_input],
                         train_true_labels,
                         validation_split=validation_split,
                         shuffle=True,
                         sample_weight=sample_weights,
                         epochs=epochs,
                         batch_size=batch_size,
                         verbose=verbose)
    
    return model, history


def evaluate_model(model, test_left_kmers_input, test_right_kmers_input, test_read_features_input, test_consensus_features_input, test_rest_features_input, test_true_labels):
    
    test_loss, test_accuracy = model.evaluate([test_left_kmers_input, test_right_kmers_input, test_read_features_input, test_consensus_features_input, test_rest_features_input],
                                              test_true_labels)
    
    return test_loss, test_accuracy


def plot_training_history(history):
    
    plt.figure()
    
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.legend(['training accuracy','validation accuracy'])
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    
    plt.savefig("training_history.png", bbox_inches='tight')
    
    return None


def generate_sample_weights(training_data, class_weights):
    
    sample_weights = [[class_weights[x] for x in y] for y in training_data]
    
    return np.array(sample_weights)


def get_keras_mask(train_rest_features_input):
    
    masking_rest_features = Masking(mask_value=-1., input_shape=train_rest_features_input.shape[1:])(train_rest_features_input)
    keras_mask = masking_rest_features._keras_mask
    
    return keras_mask


if __name__ == "__main__":

    train_left_kmers_input_file = "train_left_kmers_input.npy"
    train_right_kmers_input_file = "train_right_kmers_input.npy"
    train_read_features_input_file = "train_read_features_input.npy"
    train_consensus_features_input_file = "train_consensus_features_input.npy"
    train_rest_features_input_file = "train_rest_features_input.npy"
    train_true_labels_file = "train_true_labels.npy"
    
    train_left_kmers_input = np.load(train_left_kmers_input_file)
    train_right_kmers_input = np.load(train_right_kmers_input_file)
    train_read_features_input = np.load(train_read_features_input_file)
    train_consensus_features_input = np.load(train_consensus_features_input_file)
    train_rest_features_input = np.load(train_rest_features_input_file)
    train_true_labels = np.load(train_true_labels_file)
    print("train_left_kmers_input shape:", train_left_kmers_input.shape)
    print("train_right_kmers_input shape:", train_right_kmers_input.shape)
    print("train_read_features_input shape:", train_read_features_input.shape)
    print("train_consensus_features_input shape:", train_consensus_features_input.shape)
    print("train_rest_features_input shape:", train_rest_features_input.shape)
    print("train_true_labels shape:", train_true_labels.shape)
    
    test_left_kmers_input_file = "test_left_kmers_input.npy"
    test_right_kmers_input_file = "test_right_kmers_input.npy"
    test_read_features_input_file = "test_read_features_input.npy"
    test_consensus_features_input_file = "test_consensus_features_input.npy"
    test_rest_features_input_file = "test_rest_features_input.npy"
    test_true_labels_file = "test_true_labels.npy"
    
    test_left_kmers_input = np.load(test_left_kmers_input_file)
    test_right_kmers_input = np.load(test_right_kmers_input_file)
    test_read_features_input = np.load(test_read_features_input_file)
    test_consensus_features_input = np.load(test_consensus_features_input_file)
    test_rest_features_input = np.load(test_rest_features_input_file)
    test_true_labels = np.load(test_true_labels_file)
    print("test_left_kmers_input shape:", test_left_kmers_input.shape)
    print("test_right_kmers_input shape:", test_right_kmers_input.shape)
    print("test_read_features_input shape:", test_read_features_input.shape)
    print("test_consensus_features_input shape:", test_consensus_features_input.shape)
    print("test_rest_features_input shape:", test_rest_features_input.shape)
    print("test_true_labels shape:", test_true_labels.shape)
        
    pickle_in = open("class_weights_dict.pickle","rb")
    class_weights_dict = pickle.load(pickle_in)
    pickle_in.close()
    print("saved-computed class_weights_dict:", class_weights_dict)    
    
    # manually set class_weights_dict
    class_weights_dict = {0: 3, 1: 3, 2: 3, 3: 3, 4: 3, 5: 1.7, 6: 1.7, 7: 1.7, 8: 1.7, 9: 1}
    print("manually-set class_weights_dict:", class_weights_dict)
    
    model = build_model(train_left_kmers_input.shape[1:], train_right_kmers_input.shape[1:], train_read_features_input.shape[1:], train_consensus_features_input.shape[1:], train_rest_features_input.shape[1:])
    
    
    #epochs = 150
    epochs = 180
    batch_size = 256
    verbose = 2
    validation_split = 0.1
    
    """
    # TESTING
    train_left_kmers_input = train_left_kmers_input[:5000]
    train_right_kmers_input = train_right_kmers_input[:5000]
    train_read_features_input = train_read_features_input[:5000]
    train_consensus_features_input = train_consensus_features_input[:5000]
    train_rest_features_input = train_rest_features_input[:5000]
    train_true_labels = train_true_labels[:5000]
    
    test_left_kmers_input = test_left_kmers_input[:1500]
    test_right_kmers_input = test_right_kmers_input[:1500]
    test_read_features_input = test_read_features_input[:1500]
    test_consensus_features_input = test_consensus_features_input[:1500]
    test_rest_features_input = test_rest_features_input[:1500]
    test_true_labels = test_true_labels[:1500]
    # END TESTING
    """
    
    # generate sample weights based on class weights and apply keras mask to sample weights
    sample_weights = generate_sample_weights(train_true_labels, class_weights_dict)
    print("sample_weights shape:", sample_weights.shape)
    
    keras_mask = get_keras_mask(train_rest_features_input)
    print("keras_mask shape:", keras_mask.shape)    
    
    sample_weights *= keras_mask.numpy()
    print("sample_weights shape (after applying mask):", sample_weights.shape)
    
    
    model, history = train_model(model, train_left_kmers_input, train_right_kmers_input, train_read_features_input, train_consensus_features_input, train_rest_features_input, train_true_labels, epochs, batch_size, validation_split, verbose, sample_weights)
    
    plot_training_history(history)

    model.save('trained_model')
    
    test_loss, test_accuracy = evaluate_model(model, test_left_kmers_input, test_right_kmers_input, test_read_features_input, test_consensus_features_input, test_rest_features_input, test_true_labels)
    print("test_loss:", test_loss)
    print("test_accuracy:", test_accuracy)
    
    

