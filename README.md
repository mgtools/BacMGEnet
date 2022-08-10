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

## Dependencies
Python: python 3

NetworkX: a python package for network based analysis 

CRISPRone: a pipeline for CRISPR-Cas system annotation, included in this repository under the CRISPRone 

## Sample usage using a toy example
Please go to examples/ to see how it works using a toy example, and what the pipelines output. 

## Outputs
The main outputs are the predicted CRISPR-Cas systems (and CRISPR arrays), putative MGEs (phages and pladmids) that have traces caught in the CRISPR-Cas systems, and putative interaction network of the genomes/metagenomes and the MGEs (in GML format).

For the toy example, the results are under examples/ folder.  

crisprone/

crisprone/crispr -- CRISPR-Cas prediction results (annotations in gff file, and predicted cas genes)

crisprone/spacer_graph/module.all.ori -- CRISPR arrays

mgenet/

mgenet/protospacer_w_multihits-greedy.gff -- identified phages with protospacer information (in gff format)

mgenet/spacer2mge.gml -- spacer and MGE network (below shows a visualization of the network in Cytoscape when [this style file](https://github.com/mgtools/BacMGEnet/blob/main/misc/style.xml) is applied)

![This is an image](https://github.com/mgtools/BacMGEnet/blob/main/misc/toynetwork.png)

In this figure, the green rectanges are the genomes (hosts), yellow ovals are phages and red ovals are plasmids. 

## Sample results
### Phocaeicola vulgatus
See results/pvulgatus

### human gut MAGs
See results/gut

### wound microbiome
See results/wound
