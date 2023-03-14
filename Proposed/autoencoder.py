#importing libraries
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras import regularizers
from tensorflow.keras import backend as K

# Autoencoder architecture
def AE(data, val, encoding_dim1, encoding_dim2, seed,
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





# Denoising Autoencoder architecture
def DAE(data, val, encoding_dim1, encoding_dim2, seed,
                lambda_act, lambda_weight, epoch, bs):
    num_in = data.shape[1]
    
    # Adding noise to the train and test data
    data_noisy = data + 0.01*np.random.normal(loc=0.0, scale=1.0, size=data.shape)
    val_noisy = val + 0.01*np.random.normal(loc=0.0, scale=1.0, size=val.shape)

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
    DAE = keras.Model(inputs=input_data, outputs=decoded)

    Encoder = keras.Model(inputs=input_data, outputs=encoded)
    DAE.compile(optimizer='Adam' , loss='mse')
    # training
    #print('training the autoencoder')
    DAE.fit(data_noisy, data, epochs=epoch, batch_size=bs, validation_data=(val, val), verbose=0)
    #autoencoder.trainable = False   #freeze autoencoder weights

    return Encoder, DAE


def sparse_regularizer(activation_matrix):
    p = 0.01
    beta = 3
    p_hat = K.mean(activation_matrix) 
  
    KL_divergence = p*(K.log(p/p_hat)) + (1-p)*(K.log(1-p/1-p_hat))
    
    sum = K.sum(KL_divergence) 
   
    return beta * sum


# Autoencoder architecture
def SAE(data, val, encoding_dim1, encoding_dim2, seed,
                lambda_act, lambda_weight, epoch, bs):
    num_in = data.shape[1]

    # this is our input placeholder
    input_data = tf.keras.Input(shape=(num_in,))
    # first encoded representation of the input
    encoded = layers.Dense(encoding_dim1, activation='relu', 
                           activity_regularizer=sparse_regularizer,
                           kernel_regularizer=regularizers.l2(lambda_weight), 
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='encoder1')(input_data)
    # second encoded representation of the input
    encoded = layers.Dense(encoding_dim2, activation='relu', 
                           activity_regularizer=sparse_regularizer,
                           kernel_regularizer=regularizers.l2(lambda_weight),
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='encoder2')(encoded)
  
    # lossy reconstruction of the input
    decoded = layers.Dense(encoding_dim1, activation='relu',
                           activity_regularizer=sparse_regularizer,
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
    SAE = keras.Model(inputs=input_data, outputs=decoded)

    Encoder = keras.Model(inputs=input_data, outputs=encoded)
    SAE.compile(optimizer='Adam' , loss='mse')
    # training
    #print('training the autoencoder')
    SAE.fit(data, data, epochs=epoch, batch_size=bs, validation_data=(val, val), verbose=0)
    #autoencoder.trainable = False   #freeze autoencoder weights

    return Encoder, SAE





# Autoencoder architecture
def OrdinaryAE(data, val, encoding_dim1, encoding_dim2, seed,
                lambda_act, lambda_weight, epoch, bs):
    num_in = data.shape[1]

    # this is our input placeholder
    input_data = tf.keras.Input(shape=(num_in,))
    # first encoded representation of the input
    encoded = layers.Dense(encoding_dim1, activation='relu', 
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='encoder1')(input_data)
    # second encoded representation of the input
    encoded = layers.Dense(encoding_dim2, activation='relu', 
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='encoder2')(encoded)
  
    # lossy reconstruction of the input
    decoded = layers.Dense(encoding_dim1, activation='relu',
                           kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                           bias_initializer = tf.keras.initializers.Constant(0.1),
                           name='decoder1')(encoded)
    # final lossy reconstruction of the input
    decoded = layers.Dense(num_in, activation='linear', 
                           kernel_initializer=tf.keras.initializers.GlorotUniform(seed), 
                           bias_initializer = tf.keras.initializers.Zeros(),
                           name='decoder2')(decoded)

    # this model maps an input to its reconstruction
    OAE = keras.Model(inputs=input_data, outputs=decoded)

    Encoder = keras.Model(inputs=input_data, outputs=encoded)
    OAE.compile(optimizer='Adam' , loss='mse')
    # training
    #print('training the autoencoder')
    OAE.fit(data, data, epochs=epoch, batch_size=bs, validation_data=(val, val), verbose=0)
    #autoencoder.trainable = False   #freeze autoencoder weights

    return Encoder, OAE
