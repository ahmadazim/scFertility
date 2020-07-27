# Can scRNA-seq data give some indication of fertility in humans?

## Objective
The advent of single-cell RNA sequencing (scRNA-seq) technology has facilitated the investigation of gene expression at the single-cell level, making it ideal for studying the how gene expression levels correlate with phenotypic data. The objective of this project was to determine whether single-cell data can be used to give some indication of sperm cell viability, and thus fertility, in humans. Specifically, can fertile sperm cells be identified among less fertile sperm cells only by using gene expression data? By identifying such a population of more fertile sperm cells in the ejaculate, a putative list of marker genes for human fertility can be developed, and a fertility score can be produced based on the expression of such marker genes. 

## Methods and Results
The data in this study were acquired from scRNA-seq data published by Hermann et al. scRNA-seq data of surgical excess normal adult human testicular tissue, gathered using a droplet-based system from 10X Genomics Inc. Principal components (PCs) were used to compute unsupervised clusters. An example of such a plot based on a t-SNE dimensionality reduction technique is shown below.
![fig1](https://github.com/ahmadazim/scFertility/blob/master/fig1.png) 

Putative marker genes were identified for each cluster. Cell populations could be re-clustered for further analysis. For example, a subset containing more mature cell types (macrophage/perivascular and Peritubular/Sertoli/Leydig cells) was examined further by re-clustering.
![fig2](https://github.com/ahmadazim/scFertility/blob/master/fig2.png)
Gene ontology analysis was also performed on the marker genes of each of the six clusters shown above; the results of this analysis suggested that scRNA-seq data can in fact be used to discern between fertile sperm cells and infertile, less developed sperm cells. Further experimentation is needed to verify whether fertility marker gene score actually correlates with true fertility. 
