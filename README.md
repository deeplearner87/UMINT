# UMINT
Maitra, C., Seal, D.B., Das, V. and De, R.K., 2022. UMINT: unsupervised neural network for single cell multi-omics integration. BioRxiv, pp.2022-04.
https://www.biorxiv.org/content/10.1101/2022.04.21.489041v1

Now published at Frontiers in Molecular Biosciences
Maitra, C., Seal, D.B., Das, V. and De, R.K., Unsupervised neural network for single cell Multi-omics INTegration (UMINT): An application to health and disease. Frontiers in Molecular Biosciences, 10, p.335.
https://www.frontiersin.org/articles/10.3389/fmolb.2023.1184748/full

## Architecture of UMINT 
![UMINT](https://user-images.githubusercontent.com/113589317/232395894-fe78cfdb-d1e4-42eb-ad76-92c5987992ae.png)


## Running UMINT
To run UMINT, import `umint.py` from the `Proposed` directory and run the function `CombinedEncoder`. All the parameters are mentioned below for better understanding. One can also follow a .ipynb file from the `Proposed` directory.

### Requirements
To run umint, one needs to install `tensorflow`, `sklearn`, `scipy` and `pandas` packages. Installation codes are as follows:
+ `pip install tensorflow`
+ `pip install scikit-learn`
+ `pip install scipy`
+ `pip install pandas`

### Parameters
All input parameters are as follows: `layer_neuron`, `mid_neuron`, `seed`, `lambda_act`, `lambda_weight`, `epoch`, `bs`
+ `data`: List of input data matrices for training. [The input data matrices should be in the form of cells x features.]
+ `val`: List of data matrices for validation. [The validation data matrices also should be in the same format as input data matrices i.e., cells x features.]
+ `layer_neuron`: List of neurons in the modality encoding layer [modality wise].
+ `mid_neuron`: Dimension onto which the data is being projected.
+ `seed`: To reproduce the results set a seed value.
+ `lambda_act`: Activity regularizer parameter. 
+ `lambda_weight`: Kernel regularizer parameter.
+ `epoch`: Total number of iteration for training.
+ `bs`: Training batch size.

### Demo
#### Code to run UMINT
To run UMINT, one needs to import the script umint (within the `Proposed` directory) first. An example is provided below. Let `x1_train` [cells x features] and `x2_train` [cells x features] be two training datasets, coming from two different omics modalities, and `x1_test` [cells x features], `x2_test` [cells x features] be their respective counterparts for validation.
```
import umint
MyEncoder, MyAE = umint.CombinedEncoder(data=[x1_train, x2_train], val=[x1_test, x2_test],
                                        layer_neuron=[128, 10], mid_neuron=64, seed=98,
                                        lambda_act=0.0001, lambda_weight=0.001, epoch=25, bs=16)
```
#### Code to find the lower dimensional embedding.
Once **UMINT** is trained, to find the latent lower dimensional embedding produced by UMINT, run the code below.
```
low = MyEncoder.predict([x1, x2]) 
```
To integrate multiple modalities please change the input accordingly. The sizes of data, val and layer_neuron must match in order to run `umint.py` successfully. 

# Authors' Information
--------------------
Chayan Maitra

Machine Intelligence Unit, Indian Statistical Institute,
203 Barrackpore Trunk Road, Kolkata 700108, India

E-mail: chayanmath25_r@isical.ac.in

Dibyendu B. Seal

Tatras Data Services Pvt. Ltd., E64, Vasant Marg, Vasant Vihar, New Delhi 110057, India.

E-mail: dbsakc@caluniv.ac.in

Vivek Das

Novo Nordisk A/S,
Novo Nordisk Park 1, 2760 M ̊aløv, Denmark

E-mail: vivekdas.0687@gmail.com

Rajat K. De

Machine Intelligence Unit, Indian Statistical Institute,
203 Barrackpore Trunk Road, Kolkata 700108, India

E-mail: rajat@isical.ac.in


# Dataset Source
--------------
https://doi.org/10.5281/zenodo.7723340



# Graphical abstract
------------------
![UMINT](https://user-images.githubusercontent.com/113589317/234973253-4dacf7d3-4302-489b-8b9c-68a63203fbb3.jpeg)

