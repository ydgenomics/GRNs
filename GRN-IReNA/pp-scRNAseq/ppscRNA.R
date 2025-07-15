# Date: 250715
# IReNA /opt/software/miniconda3/envs/IReNA/bin/R
# Ref: https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing

library(IReNA)
library(Seurat)
library(optparse)

option_list <- list(
  make_option(c("-m", "--motif_txt"), type = "character",
              default = "/data/work/SCPipelines/pp_scRNA/tf_blinding_motif_ga2.txt",
              help = "Motif txt file [default: %default]"),
  make_option(c("-r", "--input_rds"), type = "character",
              default = "/data/work/SCPipelines/bulk_RNA_scRNA_singleR/split/seu_day-2.rds",
              help = "Input RDS file [default: %default]"),
  make_option(c("-n", "--n_hvg"), type = "integer", default = 3000,
              help = "Number of HVGs [default: %default]"),
  make_option(c("-M", "--method"), type = "character", default = "cluster",
              help = "Method: cluster or pseudotime [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

motif_txt <- opt$motif_txt
input_rds <- opt$input_rds
n_hvg <- opt$n_hvg
method <- opt$method

# ###call Mus musculus motif database
# motif1 <- Tranfac201803_Mm_MotifTFsF; head(motif1)
# ###call Homo sapiens motif database
# motif1 <- Tranfac201803_Hs_MotifTFsF; head(motif1)
# ###call Zebrafish motif database
# motif1 <- Tranfac201803_Zf_MotifTFsF; head(motif1)
# ###call Chicken motif database
# motif1 <- Tranfac201803_Ch_MotifTFsF; head(motif1)
motif1 <- read.table(motif_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE) # added by yd
head(motif1) # added by yd

# Step 1: Calculate the pseudotime
seurat_object <- readRDS(input_rds); seurat_object; head(rownames(seurat_object), n=10)

# Assay RNA changing from Assay5 to Assay # added by yd
seurat_object[["RNA"]] <- as(seurat_object[["RNA"]], "Assay") # added by yd

method <- "cluster"   # 可改成 "pseudotime" 或其他值做测试

if (method == "cluster") {
  # ====== 当 method 为 "cluster" 时运行的代码 ======
  cat("Running cluster analysis ...\n")
  hvg_object <- seurat_object
  hvg_object <- FindVariableFeatures(hvg_object, selection.method = "vst", nfeatures = n_hvg)
  variable_genes <- VariableFeatures(hvg_object)
  counts_matrix <- GetAssayData(hvg_object, slot = "counts")[variable_genes, ] # 提取 counts 矩阵（仅保留 variable features）
  hvg_object <- CreateSeuratObject(counts = counts_matrix)
  hvg_object@meta.data <- seurat_object@meta.data
  hvg_object <- NormalizeData(hvg_object)
  hvg_object[["RNA"]]<-as(object=hvg_object[["RNA"]],Class="Assay")
  hvg_object
  saveRDS(hvg_object, file = paste0(basename(input_rds), "_hvg", n_hvg,".rds"))
} else if (method == "pseudotime") {
  # ====== 当 method 为 "pseudotime" 时运行的代码 ======
  cat("Running pseudotime analysis ...\n")
  ### Calculate the pseudotime and return monocle object
  monocle_object <- get_pseudotime(seurat_object,gene.use = rownames(seurat_object)) #https://rdrr.io/github/jiang-junyao/IReNA/man/get_pseudotime.html
  # monocle_object <- get_pseudotime(seurat_object,gene.use = NULL)
  print(monocle_object)
  ###Add pseudotime to the Seurat object
  ### This function only support monocle object from monocle2
  seurat_with_time <- add_pseudotime(seurat_object, monocle_object)
  head(seurat_with_time@meta.data[c("Pseudotime", "State")]) # added by yd

  # Step 2: Identify DEGs and expressed transcription factors(TFs)
  ### Identify DEGs across pseudotime (qvalue < 0.05 and num_cells_expressed > 0.1)
  library(monocle)
  monocle_object <- detectGenes(monocle_object, min_expr = 1)
  monocle_object <- estimateDispersions(monocle_object)
  diff1 <- monocle::differentialGeneTest(monocle_object,fullModelFormulaStr = "~Pseudotime",relative_expr = TRUE)
  sig_genes <- subset(diff1, qval < 0.05)
  sig_genes <- subset(sig_genes, num_cells_expressed > 0.1); head(sig_genes)
  ### select Candidate TFs.
  Candidate_TFs <- c()
  for (i in 1:nrow(motif1)) {
    gene1 <- strsplit(motif1[i,5],';')[[1]]
    Candidate_TFs <- c(Candidate_TFs,gene1)
  }
  ### Identify expressed TFs
  ### Canidate TFs in our motif database are ensemble ID, If your gene names in 
  ### seurat object are Symbol ID, you need to transfer 
  ### Canidate TFs to Symbol ID first.
  expressed_tf <- rownames(extract_expressed_TFs(seurat_object,Candidate_TFs))
  expressed_tf <- expressed_tf[!expressed_tf%in%rownames(sig_genes)]
  ### Refine the seurat object
  refined_seurat <- subset(seurat_with_time, features = c(expressed_tf,rownames(sig_genes)))

  saveRDS(refined_seurat, file = paste0(basename(input_rds), "_seurat_with_time.rds"))
  # # if you already have identified DEGs, you just need to run subset function in seurat:
  # ### DEGs used here is the character class
  # seurat_with_time <- subset(seurat_with_time, features = DEGs)
} else {
  # ====== 其余情况 ======
  cat("Unknown method:", method, "\n")
}