#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# DROPSEQ Analysis to accompany:
#
# Bacon et al, 2018      
# Single-Cell Analysis Identifies Thymic Maturation Delay in Growth-Restricted Neonates 
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
 library("ggraph")
 library("cowplot")
 library("Seurat")
 library("Matrix")
 library("useful")
 library("reshape2")
 library("biomaRt")
 library("scran")
 library("scater")
 library("viridis")
 library("clustree")
 library("SingleCellExperiment")
 library("pheatmap")
 library("tidyr")
})

genecut       <- 300
mincells      <- 3
#resolution    <- "0.6"
resolution    <- "0.4_to_1.0"
slxID         <- "SLX-7632"
normalisation <- "log2"
perplexity    <- 50
methodRun     <- "m1"; 
regression    <- "umi_mt_cc"

message(paste0("Experiment      = ", slxID))
message(paste0("Gene Cut        = ", genecut))
message(paste0("Min Cells       = ", mincells))
message(paste0("Perplexity      = ", perplexity))
message(paste0("Resolution      = ", resolution)) 
message(paste0("Normalisation   = ", normalisation))
message(paste0("Method Run      = ", methodRun))
message(paste0("Regression type = ", regression))

baseDir <- "/Users/rhamilto/Documents/CTR-DROPSEQ/CTR_DropSeq_0009/SplitTests"
setwd(baseDir)

message("+--- Set up some useful functions --+")

centroidFunction <- function(tSNEClusterTable){ 
  return (data.frame( x=mean(tSNEClusterTable$tSNE_1),y=mean(tSNEClusterTable$tSNE_2))) }

meanCentroidDist <- function(tSNEClusterTable, centroid){ 
  median <- median( apply(tSNEClusterTable,1,function(x,centroid) 
                    {(sqrt((tSNEClusterTable$tSNE_1 - centroid[1])^2+(tSNEClusterTable$tSNE_2-centroid[2])^2))},centroid) ) 
  return( round(median, 3) ) }



message("+--- Read in the Robj file ---+")

if( (methodRun == "m1") & (regression == "umi_mt_cc")){
  load(paste0("SLX-7632.matrix.su.WTvar.genecut.300.res.",   resolution, ".reg.umi_mt_cc.robj"), verbose=TRUE)
  load(paste0("SLX-7632.matrix.WTvar.tsne.genecut.300.res.", resolution, ".reg.umi_mt_cc.robj"), verbose=TRUE)
  matrix.su       <- matrix.su.WTvar
  matrix.tsne     <- matrix.WTvar.tsne
  unlink("matrix.su.WTvar")
  unlink("matrix.tsneWTvar")
} else {
  message("ERROR: No matching expt details (methodRun and/or regression")
}

message("+--- Set up some colour palets for consistency ---+")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cluster_col_range    = gg_color_hue(8)
experiment_col_range = gg_color_hue(2)


message("+---------------------------------------------------------------------------------+")
message("+--- Summary of genes / cells in the matrix                                    ---+")
message("+---------------------------------------------------------------------------------+")


matrix.su@meta.data$Experiment <- gsub("\\.[0-9]$",   "",          rownames(matrix.su@meta.data))
matrix.su@meta.data$Experiment <- gsub("..{12}$","",               matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- gsub("SLX......","",             matrix.su@meta.data$Experiment)

matrix.su@meta.data$Nextera    <- matrix.su@meta.data$Experiment

matrix.su@meta.data$Experiment <- gsub("N701|N704|N705|N706","WT", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- gsub("N702|N703|N707",     "KO", matrix.su@meta.data$Experiment)

numrows  <- nrow(matrix.su@data)
numcells <- ncol(matrix.su@data)

message( paste("+---  Number of cells=", numcells, " Number of genes=", numrows, " ---+", sep=""))
if(numcells < 100){ print("+---  Below 100 cells after filtering ---+"); return(numcells)}



message("+---  nUMI Summary (ALL) ---+")
summary( matrix.su@meta.data$nUMI )
message("+---  nUMI Summary (WT) ---+")
summary( subset(matrix.su@meta.data, Experiment  == "WT")$nUMI )
message("+---  nUMI Summary (KO) ---+")
summary( subset(matrix.su@meta.data, Experiment  == "KO")$nUMI )


message("+---  nGene Summary (All) ---+")
summary( matrix.su@meta.data$nGene )
message("+---  nGene Summary (WT) ---+")
summary( subset(matrix.su@meta.data, Experiment  == "WT")$nGene )
message("+---  nGene Summary (KO) ---+")
summary( subset(matrix.su@meta.data, Experiment  == "KO")$nGene )



table(matrix.su@meta.data$Experiment, matrix.su@meta.data$Nextera)

table(matrix.su@meta.data$res.0.6, matrix.su@meta.data$Nextera )

table(matrix.su@meta.data$res.0.6, matrix.su@meta.data$Experiment )


# Reads Per Cell ALL
summary( colSums(matrix.su@raw.data[, rownames(matrix.su@meta.data) ]) )
                          
# Reads Per Cell WT
summary( colSums(matrix.su@raw.data[, rownames(subset(matrix.su@meta.data, Experiment  == "WT")) ]) )
                          
# Reads Per Cell KO
summary( colSums(matrix.su@raw.data[, rownames(subset(matrix.su@meta.data, Experiment  == "KO")) ]) )                          



files <- list.files(path=paste0(baseDir, "/UMIStats"), pattern="*star_gene_exon_tagged.gatheredBarcodeDist.gene_umi_summary.txt", full.names=T, recursive=FALSE)

message(paste0("nextera", ",", "sum.umi", ",", "mean.umi", ",",  "sum.umi.uniq", ",", "mean.umi.uniq", ",", "nrow(umi_summary)"))

lapply(files, function(x) {
  nextera               <- gsub(".*SLX.7632.*.N7", "N7", x)
  nextera               <- gsub(".star.*","", nextera)
  umi_summary           <- read.table(x, header=T) 
  rownames(umi_summary) <- gsub("^", paste0("SLX.7632.", nextera, "."), umi_summary$CellBarcodeSeq )
  umi_summary           <- umi_summary[rownames(matrix.su@meta.data), ] 
  umi_summary           <- umi_summary %>% drop_na()
  sum.umi               <- Reduce('+', umi_summary$UMI)
  mean.umi              <- sum.umi/nrow(umi_summary)
  
  sum.umi.uniq          <- Reduce('+', umi_summary$UMIuniq)
  mean.umi.uniq         <- sum.umi.uniq/nrow(umi_summary)
  message(paste0(nextera, ",", sum.umi, ",", mean.umi, ",",  sum.umi.uniq, ",", mean.umi.uniq, ",", nrow(umi_summary)))
})





message("+---------------------------------------------------------------------------------+")
message("+--- Multi-resolution Plots                                                    ---+")
message("+---------------------------------------------------------------------------------+")

for(res in c( 0.2, 0.4, 0.6, 0.8, 1.0 ))
   {
     message(paste0("Creating tsne at ", res, " resolution"))
     matrix.tsne  <- as.data.frame(GetCellEmbeddings(object=matrix.su, reduction.type="tsne", dims.use=1:2))
     matrix.tsne$Experiment <- gsub("\\.[0-9]$",   "",          rownames(matrix.tsne))
     matrix.tsne$Experiment <- gsub("..{12}$","",               matrix.tsne$Experiment)
     matrix.tsne$Experiment <- gsub("SLX......","",             matrix.tsne$Experiment)
     matrix.tsne$Experiment <- gsub("N701|N704|N705|N706","WT", matrix.tsne$Experiment)
     matrix.tsne$Experiment <- gsub("N702|N703|N707",     "KO", matrix.tsne$Experiment)
     
     reso                <- paste0("res.", res)
     clust               <- matrix.su@meta.data %>% dplyr::select(contains(reso))
     colnames(clust)     <- c("Cluster")
     matrix.tsne$Cluster <- clust$Cluster
     
     tsne.X <- ggplot(matrix.tsne, aes(x=tSNE_1, y=tSNE_2, colour=Cluster)) +
               geom_point(alpha=0.5, size=2) +
               ggtitle(paste0(slxID, " tsne res=", resolution, " genecut=", genecut, " res=", res,
                              "\n regression=", regression, " method=", methodRun, " perplexity=", perplexity))
     clusters <- levels(factor(matrix.tsne$Cluster))
     for(i in clusters){
       centroid          <- centroidFunction(matrix.tsne[matrix.tsne$Cluster == i,])
       dispersion        <- meanCentroidDist(matrix.tsne[matrix.tsne$Cluster == i,], centroid)
       cluster_num_cells <- nrow( matrix.tsne[matrix.tsne$Cluster == i,] )
       tsne.X            <- tsne.X + annotate("rect", xmin = centroid$x-7.5, xmax = centroid$x+7.5,
                                              ymin = centroid$y-5, ymax = centroid$y+5, alpha = .5,
                                              colour='black', fill='white') +
                                     annotate("text", x = centroid$x, y = centroid$y,
                                              label = paste0("C_",i,"_",cluster_num_cells,
                                                             "\n",dispersion),size=5) 
       }
     pdf(paste0(slxID, ".", "tSNEs.colbyexpt.reg.", regression, ".method.", methodRun, 
                ".perplexity.", perplexity, ".res.", res,  ".pdf") ,width=15,height=15, onefile=FALSE)
     par(bg=NA)
     print({ tsne.X })
     dev.off()
}



message("+--- clustree packaged version ---+")
ctree <- clustree(matrix.su)
pdf(paste0(slxID, ".", "clustree.", regression, ".", resolution, ".", methodRun, ".pdf") ,width=10,height=10, onefile=FALSE)
par(bg=NA)
ctree
dev.off()


for(markers in c("H19", "Cd4", "Cd8a", "Cd8b1", "Il2ra", "Cd44", "Itm2a", "Ptprc", "Bcl11b")) # "Hba-a1"))
{
  message(paste0("Creating clustree with ", markers))
  pdf(paste0(slxID, ".", "clustree.", regression, ".", resolution, ".", methodRun, ".", markers, ".pdf") ,width=10,height=10, onefile=FALSE)
  par(bg=NA)
  print({ clustree(matrix.su, node_colour = markers, node_colour_aggr = "median") })
  dev.off()
}











message("+---------------------------------------------------------------------------------+")
message("+---0.6 resolution Plots                                                       ---+")
message("+---------------------------------------------------------------------------------+")


matrix.su    <- FindClusters(object=matrix.su, reduction.type = "pca", resolution = 0.6, print.output = FALSE, force.recalc=TRUE) 
matrix.su    <- RunTSNE(matrix.su, dims.use = 1:10, do.fast = T)

matrix.tsne            <- as.data.frame(GetCellEmbeddings(object=matrix.su, reduction.type="tsne", dims.use=1:2))
matrix.tsne$Experiment <- gsub("\\.[0-9]$",   "",          rownames(matrix.tsne))
matrix.tsne$Experiment <- gsub("..{12}$","",               matrix.tsne$Experiment)
matrix.tsne$Experiment <- gsub("SLX......","",             matrix.tsne$Experiment)
matrix.tsne$Experiment <- gsub("N701|N704|N705|N706","WT", matrix.tsne$Experiment)
matrix.tsne$Experiment <- gsub("N702|N703|N707",     "KO", matrix.tsne$Experiment)

reso                    <- paste0("res.0.6")
clust                   <- matrix.su@meta.data %>% dplyr::select(contains(reso))
colnames(clust)         <- c("Cluster")
matrix.tsne$Cluster     <- clust$Cluster

cellcycle               <- matrix.su@meta.data %>% dplyr::select(contains("Cell_cycle"))
colnames(cellcycle)     <- c("Cell_cycle")
matrix.tsne$Cell_cycle  <- cellcycle$Cell_cycle

matrix.tsne$Cluster <- gsub("^",  "x", matrix.tsne$Cluster)

matrix.tsne$Cluster <- gsub("x4", "1", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x0", "2", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x1", "3", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x3", "4", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x5", "5", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x2", "6", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x7", "7", matrix.tsne$Cluster)
matrix.tsne$Cluster <- gsub("x6", "8", matrix.tsne$Cluster)

matrix.tsne$Experiment <- gsub("KO", "P0", matrix.tsne$Experiment)


matrix.tsne.csv         <- matrix.tsne
matrix.tsne.csv$Nextera <- gsub("SLX.7632.", "", rownames(matrix.tsne.csv) )
matrix.tsne.csv$Nextera <- gsub("..{12}$","",    matrix.tsne.csv$Nextera)
matrix.tsne.csv         <- within(matrix.tsne.csv, rm("tSNE_1", "tSNE_2"))
write.csv(matrix.tsne.csv, file = paste0(slxID, ".", "Table.CellsByNexteraCellcycleExperiment.", regression, ".", methodRun, ".csv"))

table(matrix.tsne.csv$Experiment,matrix.tsne.csv$Nextera)
table(matrix.tsne.csv$Cell_cycle,matrix.tsne.csv$Nextera)
table(matrix.tsne.csv$Cluster,matrix.tsne.csv$Nextera)



message("+--- Giant Heatmap ---+")

clustertableAll <- FindAllMarkers(matrix.su, only.pos = FALSE, do.print = TRUE, print.bar=TRUE)
head(clustertableAll)

# Explore the ribosomal protein content
clustertableAll.RP <- subset(clustertableAll,    grepl(glob2rx("Rp*") , gene) )
clustertableAll.RP <- subset(clustertableAll.RP, abs(avg_logFC) > 0.7) 

# Extract the top 20 differential genes per cluster
clustertableAll.top20 <- as.data.frame( clustertableAll %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_logFC)) )

write.csv(clustertableAll.top20, file = "T-Cell.Table.Supp2.csv")


clustertableAll.top20.RP <- subset(clustertableAll.top20,  grepl(glob2rx("Rp*") , gene) )


min.ave_logFC <- min(abs(clustertableAll.top20$avg_logFC))

clustertableAll.original <- clustertableAll

# Remove RBCs = cluster 6
clustertableAll <- subset(clustertableAll, cluster != 6)
# Remove Macrophages = cluster 7
clustertableAll <- subset(clustertableAll, cluster != 7)
unique(clustertableAll$cluster)

#Get the top 20 gene names per cluster
markers.use <- as.data.frame( clustertableAll %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_logFC)) )$gene
length(markers.use)
length(unique(markers.use))


message(paste0( "log2 FC range: ", min(clustertableAll$avg_logFC), " to ", max(clustertableAll$avg_logFC)))

clusterAll.data.annot           <- as.data.frame(  matrix.su@meta.data[ , c("Cell_cycle", "res.0.6")]) 
colnames(clusterAll.data.annot) <- c("Cell_cycle", "Cluster")
clusterAll.data.annot$Genotype  <- gsub("\\.[0-9]$",   "",          rownames(clusterAll.data.annot))
clusterAll.data.annot$Genotype  <- gsub("..{12}$","",               clusterAll.data.annot$Genotype)
clusterAll.data.annot$Genotype  <- gsub("SLX......","",             clusterAll.data.annot$Genotype)
clusterAll.data.annot$Genotype  <- gsub("N701|N704|N705|N706","WT", clusterAll.data.annot$Genotype)
clusterAll.data.annot$Genotype  <- gsub("N702|N703|N707",     "P0", clusterAll.data.annot$Genotype)

#
# Apply Wendi's cluster numbering
#
clusterAll.data.annot$Cluster <- gsub("^",  "x", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x4", "1", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x0", "2", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x1", "3", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x3", "4", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x5", "5", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x2", "6", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x7", "7", clusterAll.data.annot$Cluster)
clusterAll.data.annot$Cluster <- gsub("x6", "8", clusterAll.data.annot$Cluster)

cells.2.keep          <- rownames(subset(clusterAll.data.annot, clusterAll.data.annot$Cluster != 8 & clusterAll.data.annot$Cluster != 7))
clusterAll.data.annot <- subset(clusterAll.data.annot, clusterAll.data.annot$Cluster != 8 & clusterAll.data.annot$Cluster != 7)
clusterAll.data.annot <- clusterAll.data.annot[c("Cell_cycle", "Genotype", "Cluster")]
head(clusterAll.data.annot)


clusterAll.data.annot$Cell_cycle <- gsub("G2M", "G2/M", clusterAll.data.annot$Cell_cycle)
head(clusterAll.data.annot)


ScaleCols <- colorRampPalette(colors = c("white","blue"))(255)
AnnotCols <- list( Cluster=c("1"=cluster_col_range[1], "2"=cluster_col_range[2], "3"=cluster_col_range[3], "4"=cluster_col_range[4], "5"=cluster_col_range[5], "6"=cluster_col_range[6]),  
                   Genotype=c(WT="#32CD32", P0="#BF3EFF"), Cell_cycle=c("G1"="#00BA38", "S"="magenta", "G2/M"="orange"))

clusterAll.scaledata <- matrix.su@scale.data[markers.use, cells.2.keep]
corner(clusterAll.scaledata)

clusterAll.data.annot.clust.ord <- clusterAll.data.annot[ with(clusterAll.data.annot, order(Genotype, Cell_cycle, decreasing=TRUE)), ]
clusterAll.data.annot.clust.ord <- clusterAll.data.annot.clust.ord[ with(clusterAll.data.annot.clust.ord, order(Cluster, decreasing=FALSE)), ]
clusterAll.data.annot.clust.ord <- clusterAll.data.annot.clust.ord[c("Cluster", "Genotype", "Cell_cycle" )]


clusterAll.scaledata <- clusterAll.scaledata[ , rownames(clusterAll.data.annot.clust.ord)] 
corner(clusterAll.scaledata)

clusterAll.scaledata  <- clusterAll.scaledata[!duplicated(clusterAll.scaledata), ]


colgaps.indi <- c( nrow(subset(clusterAll.data.annot, Cluster==1)), nrow(subset(clusterAll.data.annot, Cluster==2)),
                   nrow(subset(clusterAll.data.annot, Cluster==3)), nrow(subset(clusterAll.data.annot, Cluster==4)), 
                   nrow(subset(clusterAll.data.annot, Cluster==5)), nrow(subset(clusterAll.data.annot, Cluster==6)) )

colgaps <- c(colgaps.indi[1], 
             colgaps.indi[1]+colgaps.indi[2], 
             colgaps.indi[1]+colgaps.indi[2]+colgaps.indi[3],
             colgaps.indi[1]+colgaps.indi[2]+colgaps.indi[3]+colgaps.indi[4], 
             colgaps.indi[1]+colgaps.indi[2]+colgaps.indi[3]+colgaps.indi[4]+colgaps.indi[5] )

message(paste0( "clusterAll.data range: ", min(clusterAll.scaledata), " to ", max(clusterAll.scaledata)))
pdf("T-Cell.Figure.3.pdf", width=27,height=25, onefile=FALSE)
par(bg=NA)
pheatmap(clusterAll.scaledata, show_colnames = FALSE, cluster_rows=TRUE, cluster_cols=FALSE, annotation=clusterAll.data.annot, annotation_colors=AnnotCols, 
         fontsize=28, fontsize_row=20, col=ScaleCols, gaps_col=colgaps, cutree_rows=8) 
dev.off()


clustertableAll <- clustertableAll.original


message("+--- TSNE WT vs KO ---+")
tsne.wtko <- ggplot(data=as.data.frame(matrix.tsne), aes(x=tSNE_1, y=tSNE_2, colour=Experiment)) +
             geom_point(alpha=0.5, size=2) + scale_colour_manual("", values = c(WT="#32CD32", P0="#BF3EFF"), limits=c("P0", "WT")) +
             geom_density_2d(data=as.data.frame(matrix.tsne), aes(x=tSNE_1, y=tSNE_2, group=Cluster), linetype='dashed', colour="black", bins=4, alpha=0.2) 

clusters <- levels(factor(matrix.tsne$Cluster))
for(i in clusters){
  centroid            <- centroidFunction(matrix.tsne[matrix.tsne$Cluster == i,])
  centroid_label_size <- 25
  
  if(i == 7){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5; centroid_label_size <- 10}
  if(i == 8){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5; centroid_label_size <- 10}
  tsne.wtko            <- tsne.wtko + annotate("text", x = centroid$x, y = centroid$y, label = paste0(i),size=centroid_label_size, alpha=0.6) 
  }
tsne.wtko            <- tsne.wtko + xlab("") + ylab("") + coord_fixed() + 
                        theme(text=element_text(size=16,  family="sans"),
                              legend.position=c(0.01, 0.2), legend.text=element_text(size=30), legend.title=element_text(size=20), line = element_blank(),
                              legend.key.height = unit(1.25, "cm"), axis.text.x=element_blank(), axis.text.y=element_blank()) + guides(colour = guide_legend(override.aes = list(size=5)))

tsne.wtko 



pdf("T-Cell.Figure.5A.pdf", width=10,height=10, onefile=FALSE)
par(bg=NA)
tsne.wtko
dev.off()



message("+--- TSNE Cell Cycle ---+")

tsne.cc <- ggplot(data=as.data.frame(matrix.tsne), aes(x=tSNE_1, y=tSNE_2, colour=Cell_cycle)) +
           geom_point(alpha=0.5, size=2) + 
           scale_colour_manual(name="", limits=c("G1", "S", "G2M"), labels=c("G1", "S", "G2/M"), values=(c("G1"="#00BA38", "S"="magenta","G2M"="orange"))) +
           geom_density_2d(data=matrix.tsne, aes(x=tSNE_1, y=tSNE_2, group=Cluster, alpha=..level..), linetype='dashed', colour="black", bins=3, alpha=0.5) 


clusters <- levels(factor(matrix.tsne$Cluster))
for(i in clusters){
  centroid            <- centroidFunction(matrix.tsne[matrix.tsne$Cluster == i,])
  centroid_label_size <- 25
  
  if(i == 7){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5; centroid_label_size <- 10}
  if(i == 8){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5; centroid_label_size <- 10}
  tsne.cc            <- tsne.cc + annotate("text", x = centroid$x, y = centroid$y, label = paste0(i),size=centroid_label_size, alpha=0.6) 
}
tsne.cc            <- tsne.cc + xlab("") + ylab("") + coord_fixed() + 
                      theme(text=element_text(size=16,  family="sans"),
                      legend.position=c(0.01, 0.2), legend.text=element_text(size=30), legend.title=element_text(size=20), line = element_blank(),
                      legend.key.height = unit(1.25, "cm"), axis.text.x=element_blank(), axis.text.y=element_blank()) + guides(colour = guide_legend(override.aes = list(size=5)))

tsne.cc 
pdf("T-Cell.Figure.4A.pdf", width=10,height=10, onefile=FALSE)
par(bg=NA)
tsne.cc
dev.off()


message("+--- TSNE Ribosomal Protein ---+")

rp.list                    <- c("Rpl37", "Rps26", "Rpl41", "Rps29", "Rps28", "Rps15a", "Rps27", "Rpl38", "Rpl36a", "Rpl23", "Rpl39", "Rps21", "Rpl22l1")
rp.cluster.ave.exp         <- AverageExpression(object = matrix.su, use.scale=TRUE)
rp.data.exp                <- as.data.frame(rp.cluster.ave.exp[rp.list,])
rp.data.exp$gene           <- rownames(rp.data.exp)
rp.data.exp.mlt            <- melt(rp.data.exp)
colnames(rp.data.exp.mlt)  <- c("gene", "cluster", "ave.exp.scale")
rp.ave.exp.per.cluster     <- aggregate(ave.exp.scale ~ cluster, rp.data.exp.mlt, mean)
RP.scaled                  <- as.data.frame(t(matrix.su@scale.data[rp.list, ]))
RP.scaled$RP.ave           <- rowMeans(RP.scaled)
matrix.rp.tsne             <- as.data.frame(matrix.tsne)
matrix.rp.tsne$RP.ave.exp  <- RP.scaled$RP.ave
matrix.rp.tsne$RP.ave.exp  <- as.numeric(as.character(matrix.rp.tsne$RP.ave.exp))


ScaleCols <- colorRampPalette(colors = c("lightgrey", "blue"))(255)

tsne.rp <- ggplot(data=matrix.rp.tsne, aes(x=tSNE_1, y=tSNE_2)) +
           geom_point(aes(colour=RP.ave.exp), alpha=0.99, size=2) +
           scale_colour_gradient2(name="Ribosomal Protein\nAve. Expression", low="white", mid=("lightsteelblue1"), high="blue", midpoint=-0.1) +
           geom_density_2d(data=matrix.rp.tsne, aes(x=tSNE_1, y=tSNE_2, group=Cluster, alpha=..level..), linetype='dashed', colour="black", bins=3, alpha=0.5) 
  
clusters <- levels(factor(matrix.rp.tsne$Cluster))
for(i in clusters){
  centroid            <- centroidFunction(matrix.rp.tsne[matrix.rp.tsne$Cluster == i,])
  centroid_label_size <- 25
  
  if(i == 7){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5; centroid_label_size <- 10}
  if(i == 8){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5; centroid_label_size <- 10}
  tsne.rp            <- tsne.rp + annotate("text", x = centroid$x, y = centroid$y, label = paste0(i),size=centroid_label_size, alpha=0.6) 
}

tsne.rp            <- tsne.rp + xlab("") + ylab("") + coord_fixed() + 
                      theme(text=element_text(size=16,  family="sans"),
                      legend.position=c(0.82,0.85), legend.text=element_text(size=15), legend.title=element_text(size=15), line = element_blank(),
                      legend.key.height = unit(1.05, "cm"), axis.text.x=element_blank(), axis.text.y=element_blank()) + 
                      guides(colour = guide_legend(reverse=T, override.aes = list( size=5)))
tsne.rp


pdf("T-Cell.Figure.4D.pdf", width=10,height=10, onefile=FALSE)
par(bg=NA)
tsne.rp 
dev.off()




message("+--- TSNE cell cycle colouring plus split by experiment ---+")

matrix.tsne.custom            <- data.frame( (matrix.tsne) ,stringsAsFactors = FALSE)
matrix.tsne.custom$Cell_cycle <- as.character(matrix.tsne.custom$Cell_cycle)
matrix.tsne.custom$wtcols     <- matrix.tsne.custom$Cell_cycle
matrix.tsne.custom$wtcols[ grepl("P0", matrix.tsne.custom$Experiment)] <- "P0" 
matrix.tsne.custom$p0cols     <- matrix.tsne.custom$Cell_cycle
matrix.tsne.custom$p0cols[ grepl("WT", matrix.tsne.custom$Experiment)] <- "WT"
head(matrix.tsne.custom)


tsne.m1.umi_mt_cc.wt.leg <- ggplot(data=as.data.frame(matrix.tsne.custom), aes(x=tSNE_1, y=tSNE_2, colour=Cell_cycle)) +
                            geom_point(alpha=0.5, size=2) + 
                            scale_colour_manual(name="Cell Cycle", limits=c("G1", "S", "G2M"), labels=c("G1", "S", "G2/M"), values=(c("G1"="#00BA38", "S"="magenta","G2M"="orange")))


tsne.m1.umi_mt_cc.wt <- ggplot(data=as.data.frame(matrix.tsne.custom), aes(x=tSNE_1, y=tSNE_2, colour=wtcols)) +
                        geom_density_2d(data=as.data.frame(matrix.tsne), aes(x=tSNE_1, y=tSNE_2, group=Cluster), linetype='dashed', colour="black", bins=10, alpha=0.2) +
                        geom_point(alpha=0.5, size=2) + coord_fixed() +
                        ggtitle("WT (With Cell Cycle Regression)") +
                        xlab("tSNE 1") + ylab("tSNE 2") + 
                        theme(legend.position="none") +
                        scale_colour_manual(name="Cell Cycle", values=(c("S"="magenta", "G1"="#00BA38","G2M"="orange", "P0"="grey"))) 
  
tsne.m1.umi_mt_cc.p0 <- ggplot(data=as.data.frame(matrix.tsne.custom), aes(x=tSNE_1, y=tSNE_2, colour=p0cols)) +
                        geom_density_2d(data=as.data.frame(matrix.tsne), aes(x=tSNE_1, y=tSNE_2, group=Cluster), linetype='dashed', colour="black", bins=10, alpha=0.2) +
                        geom_point(alpha=0.5, size=2) + coord_fixed() +
                        ggtitle("P0 (With Cell Cycle Regression)") +
                        xlab("tSNE 1") + ylab("tSNE 2") + 
                        theme(legend.position="none") +
                        scale_colour_manual(name="Cell Cycle", values=(c("S"="magenta", "G1"="#00BA38","G2M"="orange", "WT"="grey" ))) 


load(paste0("SLX-7632.matrix.su.WTvar.genecut.300.res.",   resolution, ".reg.umi_mt.robj"), verbose=TRUE)
load(paste0("SLX-7632.matrix.WTvar.tsne.genecut.300.res.", resolution, ".reg.umi_mt.robj"), verbose=TRUE)
matrix.m1.umi_mt.su   <- matrix.su.WTvar
matrix.m1.umi_mt.tsne <- matrix.WTvar.tsne
unlink("matrix.su.WTvar")
unlink("matrix.tsneWTvar")

matrix.m1.umi_mt.tsne$Experiment  <- gsub("\\.[0-9]$",   "",          rownames(matrix.m1.umi_mt.tsne))
matrix.m1.umi_mt.tsne$Experiment  <- gsub("..{12}$","",               matrix.m1.umi_mt.tsne$Experiment)
matrix.m1.umi_mt.tsne$Experiment  <- gsub("SLX......","",             matrix.m1.umi_mt.tsne$Experiment)
matrix.m1.umi_mt.tsne$Experiment  <- gsub("N701|N704|N705|N706","WT", matrix.m1.umi_mt.tsne$Experiment)
matrix.m1.umi_mt.tsne$Experiment  <- gsub("N702|N703|N707",     "KO", matrix.m1.umi_mt.tsne$Experiment)
reso                              <- paste0("res.0.6")
clust.m1.umi_mt                   <- matrix.m1.umi_mt.su@meta.data %>% dplyr::select(contains(reso))
colnames(clust.m1.umi_mt)         <- c("Cluster")
matrix.m1.umi_mt.tsne$Cluster     <- clust.m1.umi_mt$Cluster
cellcycle.m1.umi_mt               <- matrix.su@meta.data %>% dplyr::select(contains("Cell_cycle"))
colnames(cellcycle.m1.umi_mt)     <- c("Cell_cycle")
matrix.m1.umi_mt.tsne$Cell_cycle  <- cellcycle.m1.umi_mt$Cell_cycle
matrix.m1.umi_mt.tsne$Cluster     <- gsub("^",  "x", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x4", "1", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x0", "2", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x1", "3", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x3", "4", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x5", "5", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x2", "6", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x7", "7", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Cluster     <- gsub("x6", "8", matrix.m1.umi_mt.tsne$Cluster)
matrix.m1.umi_mt.tsne$Experiment  <- gsub("KO", "P0", matrix.m1.umi_mt.tsne$Experiment)
matrix.m1.umi_mt.tsne             <- data.frame( (matrix.m1.umi_mt.tsne) ,stringsAsFactors = FALSE)
matrix.m1.umi_mt.tsne$Cell_cycle  <- as.character(matrix.m1.umi_mt.tsne$Cell_cycle)
matrix.m1.umi_mt.tsne$wtcols      <- matrix.m1.umi_mt.tsne$Cell_cycle
matrix.m1.umi_mt.tsne$wtcols[ grepl("P0", matrix.m1.umi_mt.tsne$Experiment)] <- "P0" 
matrix.m1.umi_mt.tsne$p0cols      <- matrix.m1.umi_mt.tsne$Cell_cycle
matrix.m1.umi_mt.tsne$p0cols[ grepl("WT", matrix.m1.umi_mt.tsne$Experiment)] <- "WT"
head(matrix.m1.umi_mt.tsne)



tsne.m1.umi_mt.wt <- ggplot(data=as.data.frame(matrix.m1.umi_mt.tsne), aes(x=tSNE_1, y=tSNE_2, colour=wtcols)) +
                     geom_density_2d(data=as.data.frame(matrix.m1.umi_mt.tsne), aes(x=tSNE_1, y=tSNE_2, group=Cluster), linetype='dashed', colour="black", bins=10, alpha=0.2) +
                     geom_point(alpha=0.5, size=2) + coord_fixed() +
                     ggtitle("WT (Without Cell Cycle Regression)") +
                     xlab("tSNE 1") + ylab("tSNE 2") + xlim(-40,50) +
                     theme(legend.position="none") +
                     scale_colour_manual(name="Cell Cycle", values=(c("S"="magenta", "G1"="#00BA38","G2M"="orange", "P0"="grey"))) 


tsne.m1.umi_mt.p0 <- ggplot(data=as.data.frame(matrix.m1.umi_mt.tsne), aes(x=tSNE_1, y=tSNE_2, colour=p0cols)) +
                     geom_density_2d(data=as.data.frame(matrix.m1.umi_mt.tsne), aes(x=tSNE_1, y=tSNE_2, group=Cluster), linetype='dashed', colour="black", bins=10, alpha=0.2) +
                     geom_point(alpha=0.5, size=2) + coord_fixed() +
                     ggtitle("P0 (Without Cell Cycle Regression)") +
                     xlab("tSNE 1") + ylab("tSNE 2") + xlim(-40,50) +
                     theme(legend.position="none") +
                     scale_colour_manual(name="Cell Cycle", values=(c("S"="magenta", "G1"="#00BA38","G2M"="orange", "WT"="grey"))) 


tsne.before_after    <- plot_grid(tsne.m1.umi_mt_cc.wt, tsne.m1.umi_mt_cc.p0, tsne.m1.umi_mt.wt, tsne.m1.umi_mt.p0, ncol = 2, nrow = 2, 
                        labels=c("A", "B", "C", "D"), label_size = 20, align = 'v', axis = 'l' )

#"A. WT (m1.umi_mt_cc)", "B. P0 (m1.umi_mt_cc)", "C. WT (m0.umi_mt)", "D. P0 (m1.umi_mt)"

legend               <- get_legend(tsne.m1.umi_mt_cc.wt.leg)



pdf("T-Cell.Figure.Supp2.pdf", width=17.5,height=15, onefile=FALSE)
par(bg=NA)
plot_grid(tsne.before_after, legend, rel_widths = c(2,0.2))
dev.off()




message("+---------------------------------------------------------------------------------+")
message("+--- TSNE Figure 2                                                             ---+")
message("+---------------------------------------------------------------------------------+")

tsne.Fig2 <- ggplot(data=matrix.tsne, aes(x=tSNE_1, y=tSNE_2, colour=Cluster)) +
             geom_point(alpha=0.5, size=1) 

clusters <- levels(factor(matrix.tsne$Cluster))
for(i in clusters){
  centroid          <- centroidFunction(matrix.tsne[matrix.tsne$Cluster == i,])
  
  if(i == 7){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5;}
  if(i == 8){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+2.5;}
  
  dispersion        <- meanCentroidDist(matrix.tsne[matrix.tsne$Cluster == i,], centroid)
  cluster_num_cells <- nrow( matrix.tsne[matrix.tsne$Cluster == i,] )
  tsne.Fig2         <- tsne.Fig2 + coord_fixed() +
                       annotate("text", x = centroid$x, y = centroid$y, label = paste0(i),size=15, alpha=0.6) 
}
tsne.Fig2         <- tsne.Fig2 + xlab("tSNE 1") + ylab("tSNE 2") + theme(legend.position="none")
tsne.Fig2

# tsne.Fig2.annot <- ggplot(data=matrix.tsne, aes(x=tSNE_1, y=tSNE_2, colour=Cluster)) +
#                    stat_density2d(aes(fill=Cluster, alpha=.01), alpha=0.2, geom = "polygon", bins=5, linetype='dashed') +
#                    geom_point(alpha=0.95, size=1) 
# 
# clusters <- levels(factor(matrix.tsne$Cluster))
# for(i in clusters){
#   centroid          <- centroidFunction(matrix.tsne[matrix.tsne$Cluster == i,])
#   dispersion        <- meanCentroidDist(matrix.tsne[matrix.tsne$Cluster == i,], centroid)
#   cluster_num_cells <- nrow( matrix.tsne[matrix.tsne$Cluster == i,] )
#   tsne.Fig2.annot   <- tsne.Fig2.annot + 
#                        annotate("text", x = centroid$x, y = centroid$y, label = paste0(i),size=10, alpha=0.95) +
#                        theme(legend.position="none") 
# }
# tsne.Fig2.annot <- tsne.Fig2.annot + xlab("tSNE 1") + ylab("tSNE 2") + 
#                    geom_curve(aes(x=-25, y=10, xend=-15, yend=-27.5), linetype="solid", 
#                                 color="black", size=2, curvature=-2.5, 
#                    arrow = arrow(length = unit(0.3, "cm")), ncp=25) +
#                    annotate("text", x = -4, y = 0, label = paste0("Differentiation"),size=10, alpha=0.99, angle = 35) +
#                    theme(line = element_blank(), text = element_blank(), title = element_blank())
# tsne.Fig2.annot

matrix.tsne.tmp <- matrix.tsne
matrix.tsne.tmp$Cluster <- gsub("1", "DN", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("2", "DP", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("3", "DP", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("4", "DP", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("5", "DP", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("6", "TMat", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("7", "Macrophage", matrix.tsne.tmp$Cluster)
matrix.tsne.tmp$Cluster <- gsub("8", "RBC", matrix.tsne.tmp$Cluster)

tsne.Fig2.annot2 <- ggplot(data=matrix.tsne.tmp, aes(x=tSNE_1, y=tSNE_2, colour=Cluster)) +
                    stat_density2d(aes(fill=Cluster, alpha=.01), alpha=0.2, geom = "polygon", bins=5, linetype='dashed') +
                    geom_point(alpha=0.95, size=0.5) 
clusters <- levels(factor(matrix.tsne.tmp$Cluster))
for(i in clusters){
  centroid          <- centroidFunction(matrix.tsne.tmp[matrix.tsne.tmp$Cluster == i,])
  
  if(i == "Macrophage"){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+3.5;}
  if(i == "RBC"){ centroid[1] <- centroid[1]-3;  centroid[2] <- centroid[2]+3.5;}  

  dispersion        <- meanCentroidDist(matrix.tsne.tmp[matrix.tsne.tmp$Cluster == i,], centroid)
  cluster_num_cells <- nrow( matrix.tsne.tmp[matrix.tsne.tmp$Cluster == i,] )
  tsne.Fig2.annot2  <- tsne.Fig2.annot2 + 
                       annotate("text", x = centroid$x, y = centroid$y, label = paste0(i),size=7.5, alpha=0.95) +
                       theme(legend.position="none") 
}

cols_merge_clust <- c( "RBC"=cluster_col_range[8], 
                       "DP"="#EE9A00", "Macrophage"=cluster_col_range[7],"DN"=cluster_col_range[1],  "TMat"=cluster_col_range[6])

tsne.Fig2.annot2 <- tsne.Fig2.annot2 + xlab("tSNE 1") + ylab("tSNE 2") + 
                    geom_curve(aes(x=-22, y=10, xend=-10, yend=-27.5), linetype="solid", 
                    color="black", size=2, curvature=-2, 
                    arrow = arrow(length = unit(0.3, "cm")), ncp=25) +
                    scale_colour_manual( values = cols_merge_clust) + 
                    scale_fill_manual( values = cols_merge_clust) + 
                    annotate("text", x = -4, y = -9, label = 'atop(italic("Differentiation"))', parse=TRUE, size=10, alpha=0.99, angle = 0) +
                    coord_fixed() +
                    theme(line = element_blank(), text = element_blank(), title = element_blank())
tsne.Fig2.annot2



message("+--- FeaturePlot and DotPlot for marker genes ---+")

markers <- c("Il2ra", "Cd8b1", "Cd8a", "Cd4", "Ccr7", "Itm2a", "Aif1", "Hba-a1")

fplt.Fig2          <- FeaturePlot(object = matrix.su, features.plot = markers, nCol=3, min.cutoff="q9", cols.use=c("lightgrey", "blue"), pt.size=1, do.return=TRUE)
fplt.Fig2$Aif1     <- fplt.Fig2$Aif1     + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$Il2ra    <- fplt.Fig2$Il2ra    + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$Cd8b1    <- fplt.Fig2$Cd8b1    + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$Cd8a     <- fplt.Fig2$Cd8a     + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$Cd4      <- fplt.Fig2$Cd4      + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$Ccr7     <- fplt.Fig2$Ccr7     + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$Itm2a    <- fplt.Fig2$Itm2a    + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
fplt.Fig2$`Hba-a1` <- fplt.Fig2$`Hba-a1` + xlab("") + ylab("") + coord_fixed() + theme(line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())


cluster.averages           <- AverageExpression(object = matrix.su, use.scale=FALSE)
my.data                    <- as.data.frame(cluster.averages[markers,])
my.data.prop               <- as.data.frame(prop.table(as.matrix(my.data), 1))
my.data.prop$gene          <- rownames(my.data.prop)
my.data.prop.mlt           <- melt(my.data.prop)
colnames(my.data.prop.mlt) <- c("gene", "cluster", "pct.exp")
head(my.data.prop.mlt,5)

cluster.ave.exp            <- AverageExpression(object = matrix.su, use.scale=TRUE)
my.data.exp                <- as.data.frame(cluster.ave.exp[markers,])
my.data.exp$gene           <- rownames(my.data.exp)
my.data.exp.mlt            <- melt(my.data.exp)
colnames(my.data.exp.mlt)  <- c("gene", "cluster", "ave.exp.scale")
head(my.data.exp.mlt,5)

my.data.combo <- merge(my.data.prop.mlt,my.data.exp.mlt, by=c("gene","cluster") )

my.data.combo$cluster <- gsub("^",  "x", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x4", "1", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x0", "2", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x1", "3", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x3", "4", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x5", "5", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x2", "6", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x7", "7", my.data.combo$cluster)
my.data.combo$cluster <- gsub("x6", "8", my.data.combo$cluster)

min(my.data.combo$ave.exp.scale)
max(my.data.combo$ave.exp.scale)
my.data.combo$ave.exp.scale <- MinMax(data = my.data.combo$ave.exp.scale, min = -1, max = 6.5)

dplt.Fig2 <- ggplot(my.data.combo[which(my.data.combo$pct.exp>0),], aes(x=gene, y=cluster, size=pct.exp, colour=ave.exp.scale)) +  
  geom_point() + 
  scale_radius(name="Proportion of\nExpression", range = c(0, 10)) +
  scale_colour_gradient2(name="Average\nExpression\n(normalised)", low="lightgrey", mid="purple", high="blue", midpoint=3.25, breaks = c( 0, 2.5, 5.0, 6.5) ) +
  scale_x_discrete(name="Known Marker Genes", limits=markers) +
  scale_y_discrete(name="Cluster Numbers", limits=c("8", "7", "6", "5", "4", "3", "2", "1")) +
  ylab("Cluster") + 
  theme(axis.text.x = element_text(angle = -45, hjust = -0.1))
  coord_fixed() +
  theme(plot.margin = unit(c(2,2,2,2), "lines"))

dplt.Fig2  

theme_set(theme_cowplot(font_size=18))
top.row.L    <- plot_grid(tsne.Fig2,        ncol = 1, nrow = 1, labels=c("A"), label_size = 20 )
top.row.M    <- plot_grid(dplt.Fig2,        ncol = 1, nrow = 1, labels=c("B"), label_size = 20 )
top.row.R    <- plot_grid(tsne.Fig2.annot2, ncol = 1, nrow = 1, labels=c("D"), label_size = 20 )
top.row      <- plot_grid(top.row.L,top.row.M, ncol = 2, nrow = 1,  rel_heights = c(1,0.8)) 
feature.row  <- plot_grid(fplt.Fig2$Il2ra, fplt.Fig2$Cd8b1, fplt.Fig2$Cd8a, fplt.Fig2$Cd4,
                          fplt.Fig2$Ccr7, fplt.Fig2$Itm2a, fplt.Fig2$Aif1, fplt.Fig2$`Hba-a1`,
                          labels=c("C", "","","","","","",""), label_size = 20, ncol = 4, nrow = 2)
 
bottom.row   <- plot_grid(top.row.R, NULL, ncol = 2, nrow = 1, rel_heights = c(1,1) )

pdf("T-Cell.Figure.2.pdf", width=10,height=15, onefile=FALSE)
par(bg=NA)
plot_grid(top.row, feature.row, bottom.row, ncol = 1, nrow = 3, rel_heights = c(1, 1, 1) )
dev.off()


message("+--- END OF SCRIPT ---+")
