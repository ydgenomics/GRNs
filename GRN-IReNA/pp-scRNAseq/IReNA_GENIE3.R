# Date: 250710
# IReNA /opt/software/miniconda3/envs/IReNA/bin/R

library(IReNA)
library(AnnotationDbi)
library(Seurat)
library(GENIE3)
library(doParallel)
library(doRNG)
library(reshape2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
motif_txt <- args[1] # /data/work/SCPipelines/pp_scRNA/tf_blinding_motif_ga.txt
input_rds <- args[2] # /data/work/SCPipelines/bulk_RNA_scRNA_singleR/split/seu_day-2.rds
n_kcluster <- as.integer(args[3]) # less than 8
n_cores <- as.integer(args[4]) # 32
if (length(args) == 7) {
  bin_method <- args[5]   # 'pseudotime' or 'ydgenomics'
  bin_key <- args[6]      # 'seurat_clusters'
  get_SmoothByBin_PseudotimeExp_yd_R <- args[7] # path
}


# Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
###Load the test data
seurat_with_time <- readRDS(input_rds)
###Get expression profiles ordered by pseudotime
if ("Pseudotime" %in% colnames(seurat_with_time@meta.data) && 
    "State" %in% colnames(seurat_with_time@meta.data)) {
  expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 50)
} else if (bin_method == "ydgenomics") {
    source(get_SmoothByBin_PseudotimeExp_yd_R)
    expression_profile <- get_SmoothByBin_PseudotimeExp_yd(
    seurat_object=seurat_with_time,
    FC = TRUE,
    Bin_key = bin_key,
    method = 'ydgenomics',
    FcType = "Q95"
  )
} else {
  stop("Unsatisfy prerequisite, couldn't product expression_profile")
}

###Filter noise and logFC in expression profile
expression_profile_filter <- filter_expression_profile(expression_profile, FC=0.01)
###K-means clustering
clustering <- clustering_Kmeans(expression_profile_filter, K1=n_kcluster)

Kmeans_clustering_ENS <- clustering
Kmeans_clustering_ENS$Symbol <- rownames(Kmeans_clustering_ENS)
Kmeans_clustering_ENS <- Kmeans_clustering_ENS[, c("Symbol", setdiff(names(Kmeans_clustering_ENS), "Symbol"))]
Kmeans_clustering_ENS[1:5,1:5]

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3","#FF7F00", "#FFFF33", "#A65628", "#F781BF")
pdf('plot_kmeans_pheatmap.pdf', width=10, height=10)
plot_kmeans_pheatmap(clustering,ModuleColor1 = colors)
dev.off()

motif1 <- read.table(motif_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE) # added by yd
head(motif1) # added by yd

# # GENIE3[optional] to infer regulatory relationships # Test it for the most small data needed time
# # Ref: https://github.com/aertslab/GENIE3
weightMat <- GENIE3(as.matrix(seurat_with_time@assays$RNA@data),nCores = n_cores) #Error in loadNamespace(x) : there is no package called ‘doParallel’
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
save(regulatory_relationships, Kmeans_clustering_ENS, motif1, file = paste0(basename(input_rds),"_GENIE3.RData"))
# load("/data/input/Files/ResultData/Workflow/W202507020005242/regulatory_relationships.RData")
# head(regulatory_relationships)
# length(regulatory_relationships$TF)