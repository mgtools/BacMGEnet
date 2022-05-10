# BacMGEnet
Inference of bacteria-phage interactions using CRISPR

This package contains computational tools for inference of bacteria-phage interactions using CRISPR arrays. 
The tools have been extensively tested on Ubuntu linux systems. 

## The tools are oganized into two main pipelines (available under pipelines/)
### crispr_ann.py
   This pipeline is for generating CRISPR arrays given a list of genomes/metagenome
### mge_net.py
   This pipeline generates bacteria-phage network using the CRISPR arrays created by crispr_ann.py

## MGE database
To run the second pipeline (mgenet.py), one needs to donwload the MGE (including phage and plasmid sequences) database and place it under the BacMGEnet/ folder or somewhere else. 

For example, after git clone the BacMGEnet to your local machine, go under the BacMGEnet folder,

mkdir mgedb

cd mgedb

wget https://omics.informatics.indiana.edu/mg/packages/mgedb.tar.gz 

tar zxvf mgedb.tar.gz

## Sample usage using a toy example
Please go to examples/ to see how it works using a toy example, and what the pipelines output. 


## Existing results
### Phocaeicola vulgatus
See results/pvulgatus

### human gut MAGs
See results/gut

### wound microbiome
See results/wound
