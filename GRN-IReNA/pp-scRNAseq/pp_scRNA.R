# Ref: https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing
# /software/miniconda/envs/IReNA/bin/R

library(IReNA)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
motif_txt <- args[1] # /data/work/SCPipelines/pp_scRNA/tf_blinding_motif_ga.txt
input_rds <- args[2] # /data/work/SCPipelines/bulk_RNA_scRNA_singleR/split/seu_day-2.rds

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