# UMINT
Maitra, C., Seal, D.B., Das, V. and De, R.K., 2022. UMINT: unsupervised neural network for single cell multi-omics integration. BioRxiv, pp.2022-04.
https://www.biorxiv.org/content/10.1101/2022.04.21.489041v1

## Architecture of UMINT 
![UMINT](https://user-images.githubusercontent.com/113589317/232395894-fe78cfdb-d1e4-42eb-ad76-92c5987992ae.png)


## Running UMINT
To run UMINT, import `umint.py` from the `Proposed` directory and run the function `CombinedEncoder`. All the parameters are mentioned below for better understanding. One can also follow a ipynb file from the `Proposed` directory.

### Requirements
To run umint, one need to install `tensorflow`, `sklearn`, `scipy` and `pandas` packages. Installation codes are as follows:
+ `pip install tensorflow`
+ `pip install scikit-learn`
+ `pip install scipy`
+ `pip install pandas`

### Parameters
All input parameters are as follows: layer_neuron, mid_neuron, seed, lambda_act, lambda_weight, epoch, bs
+ `data`: List of input data matrices for training.
+ `val`: List of data matrices for validation.
+ `layer_neuron`: List of neurons in the modality encoding layer (modality wise)
+ `mid_neuron`: Dimension onto which the data is being projected.
+ `seed`: To reproduce the results set a seed value.
+ `lambda_act`: Activity regularizer parameter. 
+ `lambda_weight`: Kernel regularizer parameter.
+ `epoch`: Total number of iteration for training.
+ `bs`: Training batch size.

### Demo
#### Code to run UMINT
To run UMINT one need to import the script umint first. The script is in the `Proposed` directory. An example is provided below. Let `x1_train` and `x2_train` be two training datasets, comming from two different omics modalities, and `x1_test`, `x2_test` be their respective counterpart for validation.
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
To integrate multiple modalities please change the input accordingly. Size of data, val and layer_neuron must match to run `umint.py`

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
![UMINT graphical abstract](https://user-images.githubusercontent.com/113589317/232399872-46be07ba-4c88-4fd6-aea8-190ff3eb5412.png)
