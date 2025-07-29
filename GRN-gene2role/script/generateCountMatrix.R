### Date: 250729 
### gene2role /opt/software/miniconda3/envs/gene2role/bin/R


library(Seurat)
library(data.table)
library(optparse)

option_list <- list(
	make_option(c("-i", "--input_rds"), type = "character", default = "/data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/bulk_RNA_scRNA_singleR/split/seu_day-2.rds",
							help = "Path to input RDS file [default: %default]"),
	make_option(c("-k", "--cluster_key"), type = "character", default = "leiden_res_0.50",
							help = "Cluster key in metadata [default: %default]"),
	make_option(c("-v", "--cluster_value"), type = "character", default = "0,1",
							help = "Cluster values separated by comma [default: %default]"),
	make_option(c("-n", "--n_hvg"), type = "integer", default = 2000,
							help = "Number of highly variable genes [default: %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
input_rds <- opt$input_rds
cluster_key <- opt$cluster_key
cluster_value <- strsplit(opt$cluster_value, split = ",")[[1]]
n_hvg <- opt$n_hvg

set.seed(123)
# setwd('/data/work/GRN-gene2role/')
#-------------------PBMC----------------
pbmc_orig <- readRDS(input_rds)
print(pbmc_orig)
pbmc_orig[["RNA"]] <- as(pbmc_orig[["RNA"]], "Assay") # added by yd
if ("celltype" %in% colnames(pbmc_orig@meta.data)) {
	pbmc_orig$celltype0 <- pbmc_orig$celltype
}
if ("orig.ident" %in% colnames(pbmc_orig@meta.data)) {
	pbmc_orig$orig.ident0 <- pbmc_orig$orig.ident
}
pbmc_orig$celltype <- pbmc_orig@meta.data[[cluster_key]]
pbmc_orig$orig.ident <- colnames(pbmc_orig)
print("---------- Seurat object meta data ----------")
print(colnames(pbmc_orig@meta.data))

if (cluster_value[1] == "all") {
	# write count
	pbmc_count <- as.data.frame(as.matrix(pbmc_orig@assays$RNA@counts))
	pbmc <- CreateSeuratObject(pbmc_count)
	pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = n_hvg)
	# Generate count matrix for multiple GRN comparison.
	pbmc_count <- pbmc_count[VariableFeatures(pbmc),] # Only keep high variable genes
	fwrite(x = pbmc_count, file = "count.csv",row.names = TRUE) # 601 X 100
	# write metadata
	pbmc_metadata <- pbmc_orig@meta.data
	write.table(pbmc_metadata, 
			file = 'metadata.csv', 
			sep = ',',
			row.names = TRUE, 
			quote = FALSE)
} else {
	pbmc_CD4_orig <- subset(pbmc_orig, pbmc_orig@meta.data[[cluster_key]] %in% cluster_value)
	pbmc_CD4_count <- as.data.frame(as.matrix(pbmc_CD4_orig@assays$RNA@counts))

	pbmc_CD4 <- CreateSeuratObject(pbmc_CD4_count)
	pbmc_CD4 <- NormalizeData(pbmc_CD4, normalization.method = "LogNormalize", scale.factor = 10000)
	pbmc_CD4 <- FindVariableFeatures(pbmc_CD4, selection.method = "vst", nfeatures = n_hvg)
	pbmc_CD4_count <- pbmc_CD4_count[VariableFeatures(pbmc_CD4),]
	fwrite(x = pbmc_CD4_count, file = "count.csv",row.names = TRUE)
	pbmc_CD4_metadata <- pbmc_CD4_orig@meta.data
	write.table(pbmc_CD4_metadata, 
							file = 'metadata.csv', 
							sep = ',',
							row.names = TRUE, 
							quote = FALSE)
}