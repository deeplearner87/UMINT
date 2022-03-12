#importing libraries
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras import regularizers

# Autoencoder architecture
def Autoencoder(data, val, encoding_dim1, encoding_dim2, seed,
                lambda_act, lambda_weight, epoch, bs):
    num_in = data.shape[1]

    # this is our input placeholder
    input_data = tf.keras.Input(shape=(num_in,))
    # first encoded representation of the input
    encoded = layers.Dense(encoding_dim1, activation='relu', 
                           activity_regularizer=regularizers.l1(lambda_act),
                           kernel_regularizer=regularizers.l2(lambda_weight), 
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='encoder1')(input_data)
    # second encoded representation of the input
    encoded = layers.Dense(encoding_dim2, activation='relu', 
                           activity_regularizer=regularizers.l1(lambda_act),
                           kernel_regularizer=regularizers.l2(lambda_weight),
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='encoder2')(encoded)
  
    # lossy reconstruction of the input
    decoded = layers.Dense(encoding_dim1, activation='relu',
                           activity_regularizer=regularizers.l1(lambda_act),
                           kernel_regularizer=regularizers.l2(lambda_weight),
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='decoder1')(encoded)
    # final lossy reconstruction of the input
    decoded = layers.Dense(num_in, activation='linear', 
                           kernel_initializer=tf.keras.initializers.GlorotUniform(seed), 
                           bias_initializer = tf.keras.initializers.Zeros(),
                           name='decoder2')(decoded)

    # this model maps an input to its reconstruction
    autoencoder = keras.Model(inputs=input_data, outputs=decoded)

    Encoder = keras.Model(inputs=input_data, outputs=encoded)
    autoencoder.compile(optimizer='Adam' , loss='mse')
    # training
    #print('training the autoencoder')
    autoencoder.fit(data, data, epochs=epoch, batch_size=bs, validation_data=(val, val), verbose=0)
    #autoencoder.trainable = False   #freeze autoencoder weights

    return Encoder, autoencoder