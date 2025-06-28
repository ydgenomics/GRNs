# Title: 01.csv.R
# Reference: https://mp.weixin.qq.com/s/9n1ITFcC3fT8uyQGlL3Qtw
# Image: GRN-SCENIC-database--01 /opt/conda/bin/R
# Author: lili, ydgenomics
# Date: 20250525
# Description: This script is the preparation of runing pySCENIC.
# Attention: Matrix comes from Seurat[RNA@counts] and the matrix is transposed.

library(Seurat)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input_rds"), type = "character", default = "input.rds", help = "Path to the input RDS file", metavar = "character"),
    make_option(c("-o", "--output_csv"), type = "character", default = "scenic.data.csv", help = "Path to the output CSV file", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

seu <- readRDS(opt$input_rds)
seu
colnames(seu@meta.data)
head(t(as.matrix(seu@assays$RNA@counts))[1:3,1:3])
head(colnames(t(as.matrix(seu@assays$RNA@counts))))

# save csv
write.csv(t(as.matrix(seu@assays$RNA@counts)),file = opt$output_csv)