# Overview

This repository contains the source code and the trained model for De novo Error Correction for Nanopore RNA-seq Long Reads Using Deep Learning.


# Installation

Download the source code and the trained model from this repository.
 
In addition, also install the following software packages and tools: <br>
[Keras-TensorFlow](https://github.com/keras-team/keras) (v 2.8.0) <br>
[SPOA](https://github.com/rvaser/spoa) (v 4.0.7) <br>
[Pychopper](https://github.com/epi2me-labs/pychopper) (v1) <br>
[cutadapt](https://github.com/marcelm/cutadapt) (v 3.5) <br>
[isONclust](https://github.com/ksahlin/isONclust) (v 0.0.6.1) <br>

Modify the directory paths (`bin_dir`) accordingly in the code of this repository (several `.py` and `.sh` files), to point to your local directory that holds the source code of this repository and your local directories that hold the executables of the corresponding software tools.


# Model

The neural network model trained with ONT RNA-seq (cDNA) data is available under the `model/` directory (`trained_model/`) in this repository.


