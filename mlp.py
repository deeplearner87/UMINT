#importing libraries
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras import regularizers

def MLP(data, label, layer_neuron, loss, output_neuron, epoch, bs):
    num_in = data[0].shape[1]
    # this is our input placeholder
    input_data = keras.layers.Input(shape=(num_in,))
    #input_data = layers.Dropout(0.2)(input_data)
    # first encoded representation of the input
    hidden_1 = layers.Dense(layer_neuron[0], activation='relu', name='hidden1')(input_data)
    #hidden_1 = layers.Dropout(0.5)(hidden_1)
    # second encoded representation of the input
    hidden_2 = layers.Dense(layer_neuron[1], activation='relu', name='hidden2')(hidden_1)
    #hidden_2 = layers.Dropout(0.5)(hidden_2)
    
    y = layers.Dense(output_neuron, activation='sigmoid', name='predictions')(hidden_2)

    classifier = keras.Model(inputs=input_data, outputs=y)
    # Compile model
    classifier.compile(loss=loss, optimizer='adam', metrics=['accuracy'])
    # Fit the model
    classifier.fit(data[0], label[0], epochs=epoch, batch_size=bs, verbose=0)

    #print('Now making predictions')
    predictions = classifier.predict(data[1])
    predictions_ = predictions.argmax(axis=1)

    return classifier, predictions_