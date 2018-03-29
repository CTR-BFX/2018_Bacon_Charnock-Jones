# Coup of T: Finding missing immune cell subtypes in growth restricted neonates using novel single-cell sequencing analysis #



Bacon, W.A. <sup>‡,1,2,3</sup>, Hamilton, R.S. <sup>‡,2,3</sup>, Kieckbusch, J. <sup>1,2</sup>, Yu, Z. <sup>4</sup>,  Abell, C. <sup>4</sup>, Colucci, F. <sup>1,2</sup> & Charnock-Jones, D.S. <sup>§,1,2</sup>

<sup>‡</sup> Co-first authors,
<sup>§</sup> Corresponding author <br>
<sup>1</sup> Department of Obstetrics & Gynaecology,
<sup>2</sup> Centre for Trophoblast Research,
<sup>3</sup> Department of Physiology, Development, & Neuroscience, University of Cambridge, Downing Site, Cambridge, CB2 3DY,
<sup>4</sup> Department of Chemistry

## Abstract ##

To be added on paper acceptance

![T-cell_tSNE](Images/T-cell_tSNE.png?raw=true=50x)


### Data Processing ###
Raw Fastq files are converted to uBAM using PicardTools:FastqToSam (v2.9.0) and aligned with STAR (v020201) to the mouse reference genome (DropSeq Version ) via DropSeqTools (v1.12) which also provides tools to perform filtering and annotation of the alignment files.

Custom tools have been developed for assessing quality control metrics in a standardised report format, determining cell gene counts and combining digital expression matrices from different sequencing runs. All tools are freely available at https://github.com/CTR-BFX/2018_Bacon_Charnock-Jones along with scripts to recreate figures X,Y,Z.

![Processing](Images/Processing.png?raw=true=50x)

Seurat (v2.0) was used to calculate PCA, find clusters, cluster gene markers and create tSNE dimensionality reduction plots.

Normalisation, batch correction and cell cycle contributions are performed with SCRAN (v1.4.5)

Pseudo-time trajectories and diffusion maps are analyses are performed with Destiny (v2.4.5)?

Resource       | URL
-------------- | --------------
Mouse Genome   | [Link](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz)
FastQC         | [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
ClusterFlow    | [DOI](http://dx.doi.org/10.12688/f1000research.10335.2)
MultiQC        | [DOI](http://dx.doi.org/10.1093/bioinformatics/btw354)


### Scripts to reproduce paper figures ###

The provided R script assumes the script is placed in a directory containing XXXXX. The script can be run interactively in R-studio or as a batch using Rscript. Note that some of the figures in the manuscript have had some label positions moved manually to prevent overlaps. R package versions are listed int he table below

Figure    | Description | Output Filename
--------- | ----------- | ------------------------
Figure 1  | Some Plot   | Figure1.pdf


R library versions used

Library   | Version | Link
--------- | ------- | ------------------------
ggplot2   | XX      | [Link](link)


### Sample Table ###

ArrayExpress or GEO submission [Link](link)

Index  | Experiment | #Cells |
------ | ---------- | ------ |
N703   | WT         | XX     |



## References ##



## Contact ##

Contact rsh46 -at- cam.ac.uk for bioinformatics related queries
