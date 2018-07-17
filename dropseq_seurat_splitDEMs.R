#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# DROPSEQ Analysis to accompany:
#
# Bacon et al, 2018      
# Single-Cell Analysis Identifies Thymic Maturation Delay in Growth-Restricted 
# Neonatal Mice
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2018_Bacon_Charnock-Jones
#
# CTR Code: CTR_DROPSEQ_0009
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics 
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------


#
# message("+--- Uncomment if packages need installing ---+")
#
#install.packages("devtools")
#library("devtools")
#curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
#sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
#install_github("satijalab/seurat")
#install.packages("dplyr", type = "source")
#devtools::install_github("lazappi/clustree", build_vignettes = TRUE)


message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
 library("dplyr")
 library("methods")
 library("utils")
 library("ggplot2")
 library("cowplot")
 library("Seurat")
 library("Matrix")
 library("useful")
 library("reshape2")
 library("biomaRt")
 library("scran")
 library("scater")
 library("SingleCellExperiment")
})

genecut       <- 300
mincells      <- 3
resolution    <- "0.4_to_1.0"
#resolution    <- "0.6"
slxID         <- "SLX-7632"
normalisation <- "log2"

# Script has to be run individually for with and without cell cycle regression
runcellcycle <- "0"


if(runcellcycle == 1){
   regression    <- "umi_mt_cc"
}else{
   regression    <- "umi_mt"
}

message(paste0("Regression    = ", regression))
message(paste0("Resolution    = ", resolution))
message(paste0("Normalisation = ", normalisation))



#
# Basedir will need updating to the appropriate location on the local system
#
baseDir <- "/storage/CTR-Projects/CTR_DropSeq/CTR_DropSeq_0009/SLX-7632/"

message("+--- Set up some useful functions --+")

centroidFunction <- function(tSNEClusterTable){ 
  return (data.frame( x=mean(tSNEClusterTable$tSNE_1),y=mean(tSNEClusterTable$tSNE_2))) }

meanCentroidDist <- function(tSNEClusterTable, centroid){ 
  median <- median( apply(tSNEClusterTable,1,function(x,centroid) 
                    {(sqrt((tSNEClusterTable$tSNE_1 - centroid[1])^2+(tSNEClusterTable$tSNE_2-centroid[2])^2))},centroid) ) 
  return( round(median, 3) ) }

message("+--- Download the ensEMBL gene identifiers ---+")
ensembl    <-  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <-  getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart = ensembl)

message("+--- Read in the DEM file (WT) ---+")
matrix.data.WT <- read.table(gzfile(
                    paste0(baseDir, "SLX-7632.XXXXXXXXXX.WT/SLX-7632.XXXXXXXXXX.WT.dge.txt.gz")), 
                    sep="\t",header=TRUE,row.names=1, stringsAsFactors = TRUE)
dim(matrix.data.WT)

message("+--- Read in the DEM file (KO) ---+")
matrix.data.KO <- read.table(gzfile(  
                    paste0(baseDir, "SLX-7632.XXXXXXXXXX.KO/SLX-7632.XXXXXXXXXX.KO.dge.txt.gz")),   
                    sep="\t",header=TRUE,row.names=1, stringsAsFactors = TRUE)
dim(matrix.data.KO)

message("+--- Read in the DEM file (WT and KO) ---+")
matrix.data <- read.table(gzfile(
                    paste0(baseDir, "SLX-7632.XXXXXXXXXX.20180117/SLX-7632.XXXXXXXXXX.dge.txt.gz")),
                    sep="\t",header=TRUE,row.names=1, stringsAsFactors = TRUE)
dim(matrix.data)



message("+--- Set Up the seutat matrix (WT) ---+")
matrix.su.WT <- CreateSeuratObject(raw.data  = matrix.data.WT, 
                                   min.cells = mincells, 
                                   min.genes = genecut, 
                                   project   = paste(slxID, "WT genecut=", genecut, " res=", resolution, sep=""))

message("+--- Set Up the seutat matrix (KO) ---+")
matrix.su.KO <- CreateSeuratObject(raw.data  = matrix.data.KO,  
                                   min.cells = mincells,   
                                   min.genes = genecut,   
                                   project   = paste(slxID, "KO genecut=", genecut, " res=", resolution, sep=""))

message("+--- wtvarplusko Set Up the seutat object for wt var genes plus KO method ---+")
matrix.su <- CreateSeuratObject(raw.data  = matrix.data,
                                min.cells = mincells,
                                min.genes = genecut,
                                project   = paste(slxID, "WTKO genecut=", genecut, " res=", resolution, sep=""))




message("+--- SCRAN/Cyclone Cell Cycle Calculations FULL WT/KO Matrix ---+")

matrix.data                    <- matrix.su@raw.data
matrix.data$external_gene_name <- rownames(matrix.data)
matrix.data                    <- merge(matrix.data, ensEMBL2id, by="external_gene_name" )
rownames(matrix.data)          <- matrix.data$ensembl_gene_id
matrix.data                    <- matrix.data[, -grep("external_gene_name", colnames(matrix.data))]
matrix.data                    <- matrix.data[, -grep("ensembl_gene_id",    colnames(matrix.data))]
matrix.data                    <- data.matrix(matrix.data, rownames.force = TRUE)
dim(matrix.data)
numcells.raw <- ncol(matrix.su@data)
message(paste0("Number of cells (processed)=", numcells.raw))
message("Number of cells from matrix.su@raw.data")
cells               <- colnames(matrix.su@raw.data)
length(cells)       
message("Number of rows from matrix.su")
length( rownames(matrix.su@data) )


if(runcellcycle == 1){
    message("+--- Cell cycle phase assignment ---+")
    sce.annot           <- SingleCellExperiment(list(counts=matrix.data))
    mm.pairs            <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    assigned            <- cyclone(sce.annot, pairs=mm.pairs)

    message( table(assigned$phases) )

    meta_data_CC            <- data.frame("Cell_cycle"=assigned$phases)
    rownames(meta_data_CC)  <- cells
    matrix.su               <- AddMetaData(matrix.su, metadata=meta_data_CC)
    meta_data_G1            <- data.frame("Cell_cycle_score_G1"=assigned$scores$G1)
    rownames(meta_data_G1)  <- cells
    matrix.su               <- AddMetaData(matrix.su, metadata=meta_data_G1)
    meta_data_S             <- data.frame("Cell_cycle_score_S"=assigned$scores$S)
    rownames(meta_data_S)   <- cells
    matrix.su               <- AddMetaData(matrix.su, metadata=meta_data_S)
    meta_data_G2M           <- data.frame("Cell_cycle_score_G2M"=assigned$scores$G2M)
    rownames(meta_data_G2M) <- cells
    matrix.su               <- AddMetaData(matrix.su, metadata=meta_data_G2M)

    matrix.su@meta.data$Cell_Cycle_Diff <- matrix.su@meta.data$Cell_cycle_score_S - matrix.su@meta.data$Cell_cycle_score_G2M

    message( head(matrix.su@meta.data$Cell_Cycle_Diff) )
  } else {
    message("+--- SKIPPING Cell cycle phase assignment ---+") 
  }

message("+--- wtvarplusko Look at MT genes (WT)---+")
mito.genes   <- grep("^mt-", rownames(matrix.su@data), value = T)
percent.mito <- colSums(expm1(matrix.su@data[mito.genes, ]))/colSums(expm1(matrix.su@data))
matrix.su    <- AddMetaData(matrix.su, percent.mito, "percent.mito")
matrix.su    <- FilterCells(object = matrix.su,
                            subset.names = c("nGene", "percent.mito"),
                            low.thresholds = c(genecut, -Inf), high.thresholds = c(2500, 0.05))

message("+--- wtvarplusko Summary of genes / cells in the matrix (WT) ---+")
numrows  <- nrow(matrix.su@data)
numcells <- ncol(matrix.su@data)
message( paste("+--- wtvarplusko Number of cells=", numcells, " Number of genes=", numrows, " ---+", sep=""))
if(numcells < 100){ message("+--- wtvarplusko Below 100 cells after filtering ---+"); return(numcells)}
message("+--- wtvarplusko nUMI Summary (WT) ---+")
summary( matrix.su@meta.data$nUMI )
message("+--- wtvarplusko nGene Summary (WT) ---+")
summary( matrix.su@meta.data$nGene )

message("+--- wtvarplusko Normalise Data (WT) ---+")
matrix.su <- NormalizeData(object = matrix.su, normalization.method = "LogNormalize", scale.factor = 1e4)

matrix.su.WTvar <- matrix.su

message("+--- wtvarplusko FindVariableGenes (WT) ---+")

matrix.su <- FindVariableGenes(object = matrix.su, do.plot=FALSE, mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

message("+--- wtvarplusko Scale Data (WT) ---+")
if(runcellcycle == 1) {
    matrix.su       <- ScaleData(object = matrix.su, vars.to.regress = c("nUMI", "percent.mito", "Cell_Cycle_Diff"))
    matrix.su.WTvar <- ScaleData(object = matrix.su.WTvar, vars.to.regress = c("nUMI", "percent.mito", "Cell_Cycle_Diff"))
} else  {
    matrix.su       <- ScaleData(object = matrix.su, vars.to.regress = c("nUMI", "percent.mito"))
    matrix.su.WTvar <- ScaleData(object = matrix.su.WTvar, vars.to.regress = c("nUMI", "percent.mito"))
  }




if(runcellcycle == 1){
  message("+--- SCRAN/Cyclone Cell Cycle Calculations ONLY WT Matrix ---+")

  matrix.WT.data                    <- matrix.su.WT@raw.data
  matrix.WT.data$external_gene_name <- rownames(matrix.WT.data)
  matrix.WT.data                    <- merge(matrix.WT.data, ensEMBL2id, by="external_gene_name" )
  rownames(matrix.WT.data)          <- matrix.WT.data$ensembl_gene_id
  matrix.WT.data                    <- matrix.WT.data[, -grep("external_gene_name", colnames(matrix.WT.data))]
  matrix.WT.data                    <- matrix.WT.data[, -grep("ensembl_gene_id",    colnames(matrix.WT.data))]
  matrix.WT.data                    <- data.matrix(matrix.WT.data, rownames.force = TRUE)
  dim(matrix.WT.data)
  numcells.raw <- ncol(matrix.su.WT@data)

  message( paste0("Number of cells (processed)=", numcells.raw))
  message("Number of cells from matrix.su.WT@raw.data")
  cells               <- colnames(matrix.su.WT@raw.data)
  length(cells)
  message("Number of rows from matrix.su.WT")
  length( rownames(matrix.su.WT@data) )
  message("+--- Cell cycle phase assignment ---+")
  sce.annot           <- SingleCellExperiment(list(counts=matrix.WT.data))
  mm.pairs            <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  assigned            <- cyclone(sce.annot, pairs=mm.pairs)

  message( table(assigned$phases) )

  meta_data_CC            <- data.frame("Cell_cycle"=assigned$phases)
  rownames(meta_data_CC)  <- cells
  matrix.su.WT               <- AddMetaData(matrix.su.WT, metadata=meta_data_CC)
  meta_data_G1            <- data.frame("Cell_cycle_score_G1"=assigned$scores$G1)
  rownames(meta_data_G1)  <- cells
  matrix.su.WT               <- AddMetaData(matrix.su.WT, metadata=meta_data_G1)
  meta_data_S             <- data.frame("Cell_cycle_score_S"=assigned$scores$S)
  rownames(meta_data_S)   <- cells
  matrix.su.WT               <- AddMetaData(matrix.su.WT, metadata=meta_data_S)
  meta_data_G2M           <- data.frame("Cell_cycle_score_G2M"=assigned$scores$G2M)
  rownames(meta_data_G2M) <- cells
  matrix.su.WT               <- AddMetaData(matrix.su.WT, metadata=meta_data_G2M)

  matrix.su.WT@meta.data$Cell_Cycle_Diff <- matrix.su.WT@meta.data$Cell_cycle_score_S - matrix.su.WT@meta.data$Cell_cycle_score_G2M
} else {
  message("+--- SKIPPING: SCRAN/Cyclone Cell Cycle Calculations ONLY WT Matrix ---+")
}

if(runcellcycle == 1){
  message("+--- SCRAN/Cyclone Cell Cycle Calculations ONLY KO Matrix ---+")

  matrix.KO.data                    <- matrix.su.KO@raw.data
  matrix.KO.data$external_gene_name <- rownames(matrix.KO.data)
  matrix.KO.data                    <- merge(matrix.KO.data, ensEMBL2id, by="external_gene_name" )
  rownames(matrix.KO.data)          <- matrix.KO.data$ensembl_gene_id
  matrix.KO.data                    <- matrix.KO.data[, -grep("external_gene_name", colnames(matrix.KO.data))]
  matrix.KO.data                    <- matrix.KO.data[, -grep("ensembl_gene_id",    colnames(matrix.KO.data))]
  matrix.KO.data                    <- data.matrix(matrix.KO.data, rownames.force = TRUE)
  dim(matrix.KO.data)
  numcells.raw <- ncol(matrix.su.KO@data)

  message( paste0("Number of cells (processed)=", numcells.raw))
  message("Number of cells from matrix.su.KO@raw.data")
  cells               <- colnames(matrix.su.KO@raw.data)
  length(cells)
  message("Number of rows from matrix.su.KO")
  length( rownames(matrix.su.KO@data) )
  message("+--- Cell cycle phase assignment ---+")
  sce.annot           <- SingleCellExperiment(list(counts=matrix.KO.data))
  mm.pairs            <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  assigned            <- cyclone(sce.annot, pairs=mm.pairs)

  message( table(assigned$phases) )

  meta_data_CC            <- data.frame("Cell_cycle"=assigned$phases)
  rownames(meta_data_CC)  <- cells
  matrix.su.KO               <- AddMetaData(matrix.su.KO, metadata=meta_data_CC)
  meta_data_G1            <- data.frame("Cell_cycle_score_G1"=assigned$scores$G1)
  rownames(meta_data_G1)  <- cells
  matrix.su.KO               <- AddMetaData(matrix.su.KO, metadata=meta_data_G1)
  meta_data_S             <- data.frame("Cell_cycle_score_S"=assigned$scores$S)
  rownames(meta_data_S)   <- cells
  matrix.su.KO               <- AddMetaData(matrix.su.KO, metadata=meta_data_S)
  meta_data_G2M           <- data.frame("Cell_cycle_score_G2M"=assigned$scores$G2M)
  rownames(meta_data_G2M) <- cells
  matrix.su.KO               <- AddMetaData(matrix.su.KO, metadata=meta_data_G2M)

  matrix.su.KO@meta.data$Cell_Cycle_Diff <- matrix.su.KO@meta.data$Cell_cycle_score_S - matrix.su.KO@meta.data$Cell_cycle_score_G2M
}else{
  message("+--- SKIPPING SCRAN/Cyclone Cell Cycle Calculations ONLY KO Matrix ---+")
}




message("+--- Look at MT genes (WT)---+")
mito.genes   <- grep("^mt-", rownames(matrix.su.WT@data), value = T)
head(mito.genes)
head( colSums(expm1(matrix.su.WT@data[mito.genes, ])) )
percent.mito <- colSums(expm1(matrix.su.WT@data[mito.genes, ]))/colSums(expm1(matrix.su.WT@data))
head(percent.mito)
matrix.su.WT <- AddMetaData(matrix.su.WT, percent.mito, "percent.mito")
matrix.su.WT <- FilterCells(object = matrix.su.WT, 
                            subset.names = c("nGene", "percent.mito"), 
                            low.thresholds = c(genecut, -Inf), high.thresholds = c(2500, 0.05))

message("+--- Look at MT genes (KO)---+")
mito.genes   <- grep("^mt-", rownames(matrix.su.KO@data), value = T)
percent.mito <- colSums(expm1(matrix.su.KO@data[mito.genes, ]))/colSums(expm1(matrix.su.KO@data))
matrix.su.KO <- AddMetaData(matrix.su.KO, percent.mito, "percent.mito")
matrix.su.KO <- FilterCells(object = matrix.su.KO,   
                            subset.names = c("nGene", "percent.mito"),   
                            low.thresholds = c(genecut, -Inf), high.thresholds = c(2500, 0.05))


message("+--- Summary of genes / cells in the matrix (WT) ---+")
numrows  <- nrow(matrix.su.WT@data)
numcells <- ncol(matrix.su.WT@data)
message( paste("+--- Number of cells=", numcells, " Number of genes=", numrows, " ---+", sep=""))
if(numcells < 100){ message("+--- Below 100 cells after filtering ---+"); return(numcells)}
message("+--- nUMI Summary (WT) ---+")
summary( matrix.su.WT@meta.data$nUMI ) 
message("+--- nGene Summary (WT) ---+")
summary( matrix.su.WT@meta.data$nGene ) 

message("+--- Summary of genes / cells in the matrix (KO) ---+")
numrows  <- nrow(matrix.su.KO@data)
numcells <- ncol(matrix.su.KO@data)
message( paste("+--- Number of cells=", numcells, " Number of genes=", numrows, " ---+", sep=""))
if(numcells < 100){ message("+--- Below 100 cells after filtering ---+"); return(numcells)}
message("+--- nUMI Summary (KO) ---+")
summary( matrix.su.KO@meta.data$nUMI )
message("+--- nGene Summary (KO) ---+")
summary( matrix.su.KO@meta.data$nGene )



message("+--- Normalise Data (WT) ---+")
matrix.su.WT <- NormalizeData(object = matrix.su.WT, normalization.method = "LogNormalize", scale.factor = 1e4)
message("+--- Scale Data (WT) ---+")
if(runcellcycle == 1){
  matrix.su.WT <- ScaleData(object = matrix.su.WT, vars.to.regress = c("nUMI", "percent.mito", "Cell_Cycle_Diff"))
}else{
  matrix.su.WT <- ScaleData(object = matrix.su.WT, vars.to.regress = c("nUMI", "percent.mito"))
}

message("+--- Find Variable genes (WT) ---+")
matrix.su.WT <- FindVariableGenes(object = matrix.su.WT, do.plot=FALSE,
                                  mean.function = ExpMean,
                                  dispersion.function = LogVMR,
                                  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)


message("+--- Normalise Data (KO) ---+")
matrix.su.KO <- NormalizeData(object = matrix.su.KO, normalization.method = "LogNormalize", scale.factor = 1e4)
message("+--- Find Variable genes (KO) ---+")
if(runcellcycle == 1){
  matrix.su.KO <- ScaleData(object = matrix.su.KO, vars.to.regress = c("nUMI", "percent.mito", "Cell_Cycle_Diff"))
}else{
  matrix.su.KO <- ScaleData(object = matrix.su.KO, vars.to.regress = c("nUMI", "percent.mito"))
}

message("+--- Find Variable genes (KO) ---+")
matrix.su.KO <- FindVariableGenes(object = matrix.su.KO, do.plot=FALSE,
                                  mean.function = ExpMean,
                                  dispersion.function = LogVMR,
                                  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)


message("+--- Find union of variable genes ---+")
hvg.matrix.su.WT <- rownames(x = head(x = matrix.su.WT@hvg.info, n = 1000))
hvg.matrix.su.KO <- rownames(x = head(x = matrix.su.KO@hvg.info, n = 1000))
hvg.union        <- union(x = hvg.matrix.su.WT, y = hvg.matrix.su.KO)

head( hvg.union)
message( length(hvg.union) )

matrix.su.WT@meta.data[, "protocol"] <- "WT"
matrix.su.KO@meta.data[, "protocol"] <- "KO"





message("+--- Method 1 : PCA with WT variable genes ---+")

message("+--- Method 1 : Run PCA ---+")
matrix.su.WTvar <- RunPCA(object = matrix.su.WTvar, pc.genes = hvg.matrix.su.WT,
                          do.print = TRUE, pcs.print = 1:5, genes.print = 5)
matrix.su.WTvar <- ProjectPCA(matrix.su.WTvar)

message("+--- Method 1 : Find Clusters ---+")
matrix.su.WTvar <- FindClusters(object=matrix.su.WTvar, dims.use = 1:10, reduction.type = "pca", 
                                resolution = 0, print.output = FALSE, save.SNN = TRUE)
for(res in c( 0.2, 0.4, 0.6, 0.8, 1.0 ))
   {
     message(paste0("+--- Method 1 : Finding Clusters at ", res, " resolution ---+"))
     matrix.su.WTvar <- FindClusters(object=matrix.su.WTvar, reduction.type = "pca", resolution = res, print.output = FALSE) 
   }

message("+--- Method 1 : tSNE ---+")
matrix.su.WTvar <- RunTSNE(matrix.su.WTvar, dims.use = 1:10, do.fast = T)

matrix.WTvar.tsne  <- as.data.frame( GetCellEmbeddings(object = matrix.su.WTvar, 
                                                                reduction.type = "tsne", dims.use = 1:2) )
matrix.WTvar.tsne$Experiment <- gsub("\\.[0-9]$",   "",          rownames(matrix.WTvar.tsne))
matrix.WTvar.tsne$Experiment <- gsub("..{12}$","",               matrix.WTvar.tsne$Experiment)
matrix.WTvar.tsne$Experiment <- gsub("SLX......","",             matrix.WTvar.tsne$Experiment)
matrix.WTvar.tsne$Experiment <- gsub("N701|N704|N705|N706","WT", matrix.WTvar.tsne$Experiment)
matrix.WTvar.tsne$Experiment <- gsub("N702|N703|N707",     "KO", matrix.WTvar.tsne$Experiment)
matrix.WTvar.tsne$Cluster    <- matrix.su.WTvar@meta.data$res.0.6


message("+--- Method 1 : Save ROBJs ---+")
save(matrix.su.WTvar,   file=paste0(slxID, ".matrix.su.WTvar.genecut.",   genecut, ".res.", resolution, ".reg.", regression, ".robj"))
save(matrix.WTvar.tsne, file=paste0(slxID, ".matrix.WTvar.tsne.genecut.", genecut, ".res.", resolution, ".reg.", regression, ".robj"))


pdf(paste0(slxID, ".", "matrix.su.WTvar.tSNEs.colbyexpt.", regression, ".pdf") ,width=15,height=15, onefile=FALSE)
par(bg=NA)
ggplot(data=as.data.frame(matrix.WTvar.tsne), aes(x=tSNE_1, y=tSNE_2, colour=Experiment)) +          
       geom_point(alpha=0.5, size=2) +   
       ggtitle(paste0(slxID, " matrix.WTvar.tsne res=", resolution, " genecut=", genecut, "regression=", regression)) 
dev.off()

pdf(paste0(slxID, ".", "matrix.su.WTvar.tSNEs.colbyclust.", regression, ".pdf") ,width=15,height=15, onefile=FALSE)
par(bg=NA)
ggplot(data=as.data.frame(matrix.WTvar.tsne), aes(x=tSNE_1, y=tSNE_2, colour=Cluster)) +
       geom_point(alpha=0.5, size=2) +
       ggtitle(paste0(slxID, " matrix.WTvar.tsne res=", resolution, " genecut=", genecut, "regression=", regression))
dev.off()


message("+--- END OF SCRIPT ---+")
