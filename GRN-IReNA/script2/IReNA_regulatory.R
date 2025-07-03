# IReNA /opt/software/miniconda3/envs/IReNA/bin/R
args <- commandArgs(trailingOnly = TRUE)
motif_txt <- args[1] # /data/work/SCPipelines/pp_scRNA/tf_blinding_motif_ga.txt
seurat_with_time_rds <- args[2] # /data/work/SCPipelines/bulk_RNA_scRNA_singleR/split/seu_day-2.rds

library(IReNA)
library(AnnotationDbi)

# Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
###Load the test data
seurat_with_time <- readRDS(seurat_with_time_rds)
###Get expression profiles ordered by pseudotime
expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 50)
###Filter noise and logFC in expression profile
expression_profile_filter <- filter_expression_profile(expression_profile, FC=0.01)
###K-means clustering
clustering <- clustering_Kmeans(expression_profile_filter, K1=4)

clustering[1:5,1:5]

Kmeans_clustering_ENS <- clustering
Kmeans_clustering_ENS$Symbol <- rownames(Kmeans_clustering_ENS)
Kmeans_clustering_ENS <- Kmeans_clustering_ENS[, c("Symbol", setdiff(names(Kmeans_clustering_ENS), "Symbol"))]
Kmeans_clustering_ENS[1:5,1:5]

pdf('plot_kmeans_pheatmap.pdf', width=10, height=10)
plot_kmeans_pheatmap(clustering,ModuleColor1 = c('#67C7C1','#5BA6DA','#FFBF0F','#C067A9'))
dev.off()

motif_txt <- motif_txt
motif1 <- read.table(motif_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE) # added by yd
head(motif1) # added by yd

# # GENIE3[optional] to infer regulatory relationships # Test it for the most small data needed time
# # Ref: https://github.com/aertslab/GENIE3
library(Seurat)
library(GENIE3)
library(doParallel)
library(doRNG)
library(reshape2)
library(dplyr)
weightMat <- GENIE3(as.matrix(seurat_with_time@assays$RNA@data),nCores = 16) #Error in loadNamespace(x) : there is no package called ‘doParallel’
weightMat <- getLinkList(weightMat)
regulation <- weightMat[weightMat[,3]>0.0002,]
### add regulation type for each gene pair
regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS,regulation)
### check whether source genes are transcription factors
motifTF <- c()
for (i in 1:nrow(motif1)) {
  TF <- strsplit(motif1[i,5],';')[[1]]
  motifTF <- c(motifTF,TF)
}
regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]
head(regulatory_relationships)
save(regulatory_relationships, file = "regulatory_relationships.RData")