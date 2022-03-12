#importing libraries
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras import regularizers

# Proposed UMINT architecture
def CombinedEncoder(data, val, layer_neuron, mid_neuron, seed, 
                    lambda_act, lambda_weight, epoch, bs):
    num_in_neuron = [i.shape[1] for i in data]
    
    inputs = [tf.keras.Input(shape=(num_in_neuron[i],)) for i in range(len(data))]

    layer1 = [layers.Dense(layer_neuron[i], activation="relu",
                      activity_regularizer=regularizers.l1(lambda_act),
                      kernel_regularizer=regularizers.l2(lambda_weight), 
                      kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                      bias_initializer = tf.keras.initializers.Constant(0.1)
                      )(inputs[i]) for i in range(len(data))]

    mid = layers.Dense(mid_neuron, activation="relu", 
                       activity_regularizer=regularizers.l1(lambda_act),
                       kernel_regularizer=regularizers.l2(lambda_weight),
                       kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                       bias_initializer = tf.keras.initializers.Constant(0.1)
                       )(tf.keras.layers.Concatenate(axis=1)(layer1))
    
    layer3 = [layers.Dense(layer_neuron[i], activation="relu",
                      activity_regularizer=regularizers.l1(lambda_act),
                      kernel_regularizer=regularizers.l2(lambda_weight), 
                      kernel_initializer=tf.keras.initializers.he_uniform(seed), 
                      bias_initializer = tf.keras.initializers.Constant(0.1)
                      )(mid) for i in range(len(data))]
    
    outputs = [layers.Dense(num_in_neuron[i], activation="linear",
                      kernel_initializer=tf.keras.initializers.GlorotUniform(seed), 
                      bias_initializer = tf.keras.initializers.Zeros()
                      )(layer3[i]) for i in range(len(data))]
    
    MyAE = keras.Model(inputs=inputs, outputs=outputs)
    MyEncoder = keras.Model(inputs=inputs, outputs=mid)

    MyAE.compile(
        loss=keras.losses.mean_squared_error,
        optimizer=keras.optimizers.Adam()
    )

    MyAE.fit(data, data,  epochs=epoch, batch_size=bs, validation_data=(val, val), verbose=0)

    return MyEncoder, MyAE