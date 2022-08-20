# get_confusion_matrix.py
#
# Name: Laura Tung
#
# Usage: python get_confusion_matrix.py <model_file>
#
# Note: Use integer true labels.
#

import sys
import numpy as np

from tensorflow import keras
from keras.models import load_model
from keras.layers import Masking

from sklearn.metrics import confusion_matrix, classification_report


class_names = [
        'Insertion',
        'Substituted_A',
        'Substituted_C',
        'Substituted_G',
        'Substituted_T',
        'Deleted_A',
        'Deleted_C',
        'Deleted_G',
        'Deleted_T',
        'No_correction'
    ]


def get_accuracy(test_true_labels, predicted_labels, keras_mask):
    
    m = keras.metrics.SparseCategoricalAccuracy()
    m.update_state(test_true_labels, predicted_labels, sample_weight=keras_mask)
    accuracy = m.result().numpy()
    
    return accuracy


def get_accuracy_simple(unmasked_test_true_labels, unmasked_integer_predicted_labels):
    
    m = keras.metrics.Accuracy()
    m.update_state(unmasked_test_true_labels, unmasked_integer_predicted_labels)
    accuracy = m.result().numpy()
    
    return accuracy


def apply_mask_into_1d(label_array, keras_mask):
    
    mask = keras_mask.numpy()
    
    x = np.ma.array(label_array, mask=~mask)
    unmasked_1d_label_array = x.compressed()
    
    return unmasked_1d_label_array


def get_unmasked_labels(test_true_labels, integer_predicted_labels, keras_mask):
    
    unmasked_test_true_labels = apply_mask_into_1d(test_true_labels, keras_mask)
    unmasked_integer_predicted_labels = apply_mask_into_1d(integer_predicted_labels, keras_mask)
    
    return unmasked_test_true_labels, unmasked_integer_predicted_labels


def get_keras_mask(test_rest_features_input):
    
    masking_rest_features = Masking(mask_value=-1., input_shape=test_rest_features_input.shape[1:])(test_rest_features_input)
    keras_mask = masking_rest_features._keras_mask
    
    return keras_mask
    
    
if __name__ == "__main__":
    
    model_file = sys.argv[1]
    
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
    
    """
    # TESTING
    test_left_kmers_input = test_left_kmers_input[:1500]
    test_right_kmers_input = test_right_kmers_input[:1500]
    test_read_features_input = test_read_features_input[:1500]
    test_consensus_features_input = test_consensus_features_input[:1500]
    test_rest_features_input = test_rest_features_input[:1500]
    test_true_labels = test_true_labels[:1500]
    # END TESTING
    """
    
    # load the model
    model = load_model(model_file)
    
    # get predictions
    predicted_labels = model.predict([test_left_kmers_input, test_right_kmers_input, test_read_features_input, test_consensus_features_input, test_rest_features_input])
    print("predicted_labels shape:", predicted_labels.shape)
    
    # model.predict() returns a numpy array that does not carry ._keras_mask anymore.
    # so we have to manually apply the Masking layer to an input in order to obtain the ._keras_mask
    keras_mask = get_keras_mask(test_rest_features_input)
    print("keras_mask shape:", keras_mask.shape)
    
    # confirm with the 'accuracy' by model.evaluate() previously
    test_accuracy = get_accuracy(test_true_labels, predicted_labels, keras_mask)
    print("test_accuracy:", test_accuracy)
    
    # convert probability predictions to integer class predictions
    integer_predicted_labels = np.argmax(predicted_labels, axis=2)
    print("integer_predicted_labels shape:", integer_predicted_labels.shape)
    
    # get unmasked true labels and unmasked integer predicted labels as 1D arrays
    unmasked_test_true_labels, unmasked_integer_predicted_labels = get_unmasked_labels(test_true_labels, integer_predicted_labels, keras_mask)
    
    # further verify the accuracy by using unmasked 1D true labels and unmasked 1D integer predicted labels
    test_accuracy_simple = get_accuracy_simple(unmasked_test_true_labels, unmasked_integer_predicted_labels)
    print("test_accuracy (simple direct way):", test_accuracy_simple)
    
    # now get the confusion matrix
    cf_matrix = confusion_matrix(unmasked_test_true_labels, unmasked_integer_predicted_labels)
    print("Confusion Matrix:")
    print(cf_matrix)
    
    # also get the classification report
    print("Classification report:")
    print(classification_report(unmasked_test_true_labels, unmasked_integer_predicted_labels, target_names=class_names))
    
    
    
