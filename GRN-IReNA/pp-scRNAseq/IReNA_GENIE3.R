# Date: 250715
# IReNA /opt/software/miniconda3/envs/IReNA/bin/R

library(IReNA)
library(AnnotationDbi)
library(Seurat)
library(GENIE3)
library(doParallel)
library(doRNG)
library(reshape2)
library(dplyr)
library(optparse)

option_list <- list(
  make_option(c("-m", "--motif_txt"), type = "character", default = "/data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_scRNA/tf_blinding_motif_ga2.txt", help = "Motif txt file path"),
  make_option(c("-i", "--input_rds"), type = "character", default = "/data/users/yangdong/yangdong_8632f88957bb4c4f85daf59edaf6b059/online/IReNA/P1/seu_day1.rds_hvg3000.rds", help = "Input RDS file path"),
  make_option(c("-k", "--n_kcluster"), type = "integer", default = 4, help = "Number of k-means clusters (less than 8)"),
  make_option(c("-c", "--n_cores"), type = "integer", default = 32, help = "Number of cores"),
  make_option(c("-b", "--bin_method"), type = "character", default = "cluster", help = "Bin method: 'pseudotime' or 'cluster'"),
  make_option(c("-K", "--bin_key"), type = "character", default = "seurat_clusters", help = "Bin key"),
  make_option(c("-s", "--get_SmoothByBin_PseudotimeExp_yd_R"), type = "character", default = "/data/users/yangdong/yangdong_8632f88957bb4c4f85daf59edaf6b059/online/IReNA/get_SmoothByBin_PseudotimeExp_yd.R", help = "Path to get_SmoothByBin_PseudotimeExp_yd R script")
)

opt <- parse_args(OptionParser(option_list = option_list))

motif_txt <- opt$motif_txt
input_rds <- opt$input_rds
n_kcluster <- opt$n_kcluster
n_cores <- opt$n_cores
bin_method <- opt$bin_method
bin_key <- opt$bin_key
get_SmoothByBin_PseudotimeExp_yd_R <- opt$get_SmoothByBin_PseudotimeExp_yd_R


# Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
###Load the test data
seurat_with_time <- readRDS(input_rds)
###Get expression profiles ordered by pseudotime
if (bin_method == "pseudotime" && "Pseudotime" %in% colnames(seurat_with_time@meta.data) && 
    "State" %in% colnames(seurat_with_time@meta.data)) {
  expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 50)
} else if (bin_method == "cluster") {
    source(get_SmoothByBin_PseudotimeExp_yd_R)
    expression_profile <- get_SmoothByBin_PseudotimeExp_yd(
    seurat_object=seurat_with_time,
    FC = TRUE,
    Bin_key = bin_key,
    method = 'cluster',
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

colors <- c("#e41a89", "#377EB8", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#00CED1")
colors <- colors[1:n_kcluster]
pdf('plot_kmeans_pheatmap.pdf', width=12, height=10)
plot_kmeans_pheatmap(clustering,ModuleColor1 = colors, show_colnames = TRUE, Show.Module = TRUE)
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