# pCoMiC

Predicted Classification of Mutations in CFTR (pCoMiC) was created as a fourth-year project for BIOL 469, Genomics, and BIOL 465, Structural Bioinformatics, at the University of Waterloo.

The pipeline aims to classify mutations in the coding sequence (CDS) of the Cystic Fibrosis Transmembrane Regulator (CFTR) gene in humans. Mutations can be categorized into six distinct classes: I–V and benign, each manifesting different phenotypic characteristics. Additional information on said classes is available here: https://www.cff.org/research-clinical-trials/types-cftr-mutations.

Importantly, pCoMiC is not intended for any clinical or experimental use and should not be used as a definitive classifier for CFTR mutations.

# Theory

The pipeline consists of two separate methods, both aiming to classify mutations: An initial probability-based model and a more complex machine learning model. The theory described below for both models is a general summary and a more in-depth description can be seen in the paper attached in the sup_data folder.

## Probability-Based Model

The probability-based model relies on analyzing the change in key amino acid properties, namely hydropathy, relative charge and hydrogen bonding capability, between the regions surrounding the mutated residue and the control residue. Significant deviations from the control state may be identified by the pipeline as CF-causing. 
Custom probability scoring functions were constructed to determine the probability of the mutation landing in a given class based off the degree of deviation in the average amino acid properties from baseline states. 

## Machine Learning Model

A multi-head attention transformer network was construct using PyTorch and its full architecture can be viewed either in the code or in the attached paper. The network was trained on 41 known mutations and tested on 13 known mutations.
The dataset was, unfortunately, incredibly limited due to manual annotation being required and restricted time.
As a quick workaround, a semi-supervised learning rule was adopted where the model was trained with 80 additional unlabeled samples with a weight ratio of 0.9:0.1 (labeled:unlabeled) representing the weight of the loss for both datasets in the combined loss through each epoch.

# Final Notes

Since pCoMiC is not intended for any clinical use, an installation guide has not been provided though you are more than welcome to download the repository and test the pipeline on your own. 
There are many packages that will have to be downloaded and, for the machine learning model, I have included the option to either train the network yourself or use the pre-trained network in the code files (see comments in network.py). 
I am open to any suggestions regarding improvements in efficiency, accuracy or readability.
