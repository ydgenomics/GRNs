# Title: 03.plot.R
# Date: 20250527
# Coder: lili, ydgenomics
# Description: Using SCENIC(R) to plot the results of pySCENIC
# Input: 
# Output:
# Image: GRN-allSCENIC--01 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/9n1ITFcC3fT8uyQGlL3Qtw

# check hdf5 of this image and specify!
Sys.setenv(HDF5_DIR="/usr/local/hdf5-1.8.13")
hdf5_dir <- Sys.getenv("HDF5_DIR")
print(hdf5_dir)
path <- Sys.getenv("PATH")
print(path)
ld_library_path <- Sys.getenv("LD_LIBRARY_PATH")
print(ld_library_path)
dyn.load("/usr/local/hdf5-1.8.13/lib/libhdf5.so")
dyn.load("/usr/local/hdf5-1.8.13/lib/libhdf5_hl.so")
library(hdf5r)
library(SCopeLoomR)

# library
library(Seurat)
library(SCopeLoomR) 
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap) 
library(data.table)
library(scRNAseq) 
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize)
library(optparse)

option_list <- list(
    make_option(c("-l", "--aucell_loom"), type = "character", default = "/data/work/0.peanut/GRN/peanut/H2014/aucell.loom",
                            help = "Path to aucell.loom file", metavar = "character"),
    make_option(c("-r", "--input_rds"), type = "character", default = "/data/work/0.peanut/convert/H2014_dataget_Anno_rename_threelayers.cg_cgn.rds",
                            help = "Path to input Seurat RDS file", metavar = "character"),
    make_option(c("-k", "--cluster_key"), type = "character", default = "cell",
                            help = "Cluster key in Seurat metadata", metavar = "character"),
    make_option(c("-a", "--assay"), type = "character", default = "RNA",
                            help = "Assay name in Seurat object", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# load data
#loom <- open_loom('/data/work/tomato/3.pySCENIC/aucell.loom')
loom <- open_loom(opt$aucell_loom)

# get_regulons 函数返回行为转录因子，列为基因的矩阵，值代表该基因是否属于转录因子的靶基因，是的话值为 1，不是则为 0
regulons_incidMat <- SCopeLoomR::get_regulons(loom,column.attr.name="Regulons")
regulons_incidMat[,1:5]
#          AL627309.1 AP006222.2 RP11-206L10.2 RP11-206L10.9
# ACAA1(+)          0          0             0             0
# ATF1(+)           0          0             0             0
# ATF3(+)           0          0             0             0
# ATF4(+)           0          0             0             0

#top2 <- head(rownames(regulons_incidMat),2)
#print(top2)
gene_list <- rownames(regulons_incidMat)
print(gene_list)

#regulonsToGeneLists 函数返回每个转录因子调控的靶基因列表
regulons <- SCENIC::regulonsToGeneLists(regulons_incidMat) 
head(regulons)

# get_regulons_AUC 返回计算好的转录因子的活性分数矩阵，行为转录因子，列为细胞
regulonAUC <- SCopeLoomR::get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
head(regulonAUC)

# get_regulon_thresholds 返回每个转录因子的 AUC 阈值
regulonAucThresholds <- SCopeLoomR::get_regulon_thresholds(loom)
head(regulonAucThresholds)

embeddings <- SCopeLoomR::get_embeddings(loom)
embeddings
close_loom(loom)

# load seurat object
#seu <- readRDS("/data/work/tomato/1.annotation/SixTime_SCT_cellannotation.rds")
seu <- readRDS(opt$input_rds)
seu <- UpdateSeuratObject(seu)
DefaultAssay(seu) <- opt$assay
DefaultAssay(seu) <- "RNA" # Because the fisrt step extracting the matrix is RNA, the downstream analysis is also RNA!
seu
colnames(seu@meta.data)

Idents(seu) <- opt$cluster_key

# 把每个转录因子的活性合并到 metadata 里面
# precoess
sub_regulonAUC <- regulonAUC[,match(colnames(seu),colnames(regulonAUC))]
# check
identical(colnames(sub_regulonAUC),colnames(seu))
# combine
seu@meta.data=cbind(seu@meta.data,t(assay(sub_regulonAUC)))

#############plot
#setwd("/data/work/tomato/3.pySCENIC/output1")

pdf("VlnPlot_regulons.pdf", height = 8, width = 12)
VlnPlot(seu,features=gene_list,pt.size=0)
dev.off()

pdf("FeaturePlot_regulons.pdf", height = 8, width = 12)
FeaturePlot(object=seu,features=gene_list)
dev.off()

pdf("RidgePlot_regulons.pdf", height = 4, width = 12)
RidgePlot(seu,features=gene_list,ncol=2)
dev.off()

pdf("DotPlot_regulons.pdf", height = 8, width = 12)
p <- DotPlot(seu, features = gene_list) + RotatedAxis()
print(p)
dev.off()

# 热图可视化
# ==============================================================================
# AVG expression
# ==============================================================================
cellClusters <- data.frame(
    row.names=colnames(seu),
    seurat_clusters=as.character(seu@meta.data[[opt$cluster_key]])) |>
dplyr::mutate(seurat_clusters=ifelse(is.na(seurat_clusters),"unkown",seurat_clusters))
head(cellClusters)
# 这段代码的意思是创建一个名为"cellClusters"的数据框，其中的行名（row.names）是来自于"seurat.data"数据集的列名（colnames），而"seurat_clusters"列则是将"seurat.data$seurat_annotations"转换为字符型数据后的结果。
# 接下来，使用dplyr包中的mutate函数对"seurat_clusters"列进行处理。如果"seurat_clusters"列中的值为NA（缺失值），则将其替换为"unknown"（未知），否则保持不变。最终将处理结果赋值给"seurat_clusters"列。

cellsPerGroup <- split(rownames(cellClusters),cellClusters$seurat_clusters)
#head(cellsPerGroup)

# remove extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))

# scale
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center=T,scale=T))

# heatmap
pdf("Heatmap_regulons.pdf", height = 8, width = 8)
Heatmap(matrix=t(regulonActivity_byGroup_Scaled[,]))
dev.off()

# 绘制所有细胞的转录因子活性热图
ht_auc <- assay(sub_regulonAUC)
ht_auc <- ht_auc[,as.character(unlist(cellsPerGroup))]
identical(colnames(ht_auc),as.character(unlist(cellsPerGroup))) # [1] TRUE
ht_auc_scale <- t(scale(t(ht_auc),center=T,scale=T))
# x = 1
lapply(seq_along(names(cellsPerGroup)),function(x){
    rep(names(cellsPerGroup)[x],length(cellsPerGroup[[x]]))
}) |> unlist() -> celltypes

#col_anno <- columnAnnotation(celltype = celltypes)

colors20 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#999999", "#6BAED6", "#3182BD", "#08519C", "#31A354", "#74C476", "#A1D99B", "#E6550D",
  "#FDAE6B", "#E6550D", "#31A354", "#756BB1"
)
unique_celltypes <- unique(celltypes)
length(unique_celltypes)

# Build a match between cell types and colors
custom_colors <- setNames(colors20[1:length(unique_celltypes)], unique_celltypes) # nolint
custom_colors
col_anno <- columnAnnotation(celltype = celltypes, col = list(celltype = custom_colors))

# plot
pdf("Heatmap_regulons_allcells.pdf", height = 4, width = 10)
Heatmap(
    matrix = ht_auc_scale,
    top_annotation = col_anno,
    col = colorRamp2(c(-2,0,2), c("#003399", "white", "#990066")),
    cluster_columns = F,
    show_row_names = T,
    show_column_names = F)

dev.off()

# 利用 AUC 计算的阈值，绘制二分类变量热图，也是文献经常见的
# ==============================================================================
# all cells heatmap
# ==============================================================================
ht_auc <- data.frame(ht_auc)
# x = 1
purrr::map_df(seq_along(regulonAucThresholds),function(x){
    tmp <- ht_auc[regulonAucThresholds[x],]
    val <- as.numeric(names(regulonAucThresholds[x]))
    tmp <- data.frame(apply(tmp, c(1,2), function(x) ifelse(x > val,1,0)))
                            return(tmp)
}) -> binary_tfs
binary_tfs <- binary_tfs[,as.character(unlist(cellsPerGroup))]
identical(colnames(binary_tfs),as.character(unlist(cellsPerGroup)))

# Random select three rows do annotation
#mark <- sample(rownames(binary_tfs),size = 3,replace = F)
#at <- match(mark,rownames(binary_tfs))
#tfs_mark = rowAnnotation(foo = anno_mark(at = at,labels=mark))
                      
# Do annotation for all genes
#mark <- rownames(binary_tfs)
#at <- match(mark, rownames(binary_tfs))
#tfs_mark <- rowAnnotation(foo = anno_mark(at = at, labels = mark))

# plot
pdf("Heatmap_Binary_regulons_allcells.pdf", height = 4, width = 10)
Heatmap(
    matrix = binary_tfs,
    name = "Binary activity of regulon",
    cluster_columns = F,
    col = c("white", "black"),
    top_annotation = col_anno,
    show_row_names = T,
    show_column_names = F)

dev.off()

# 对细胞聚类
pdf("Heatmap_Binary_regulons_cluster_allcells.pdf", height = 4, width = 10)
Heatmap(
    matrix = binary_tfs,
    name = "Binary activity of regulon in cluster",
    cluster_columns = T, # if cell more than 65536 will error!
    col = c("white", "black"),
    top_annotation = col_anno,
    show_row_names = T,
    show_column_names = F)

dev.off()