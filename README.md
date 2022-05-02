# BacMGEnet
Inference of bacteria-phage interactions using CRISPR

This package contains computational tools for inference of bacteria-phage interactions using CRISPR arrays. 

Please go to examples to see how it works using a toy example.  

The tools are oganized into two main pipelines (available under pipelines/)
## crispr_ann.py
   This pipeline is for generating CRISPR arrays given a list of genomes/metagenome
## mge_net.py
   This pipeline generates bacteria-phage network using the CRISPR arrays created by crispr_ann.py
