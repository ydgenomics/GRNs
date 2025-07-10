# IReNA /opt/software/miniconda3/envs/IReNA/bin/R
library(IReNA)
library(AnnotationDbi)

args <- commandArgs(trailingOnly = TRUE)
time_rds <- args[1]
input_bams_txt <- args[2]
bam_inputs <- scan(input_bams_txt, what = character(), sep = ",")
tf2motif_txt <- args[3]
genie3_rdata <- args[4]
footprints_bed <- args[5]
genome_fa <- args[6]
pwm_dir <- args[7]
peaks_bed <- args[8]
txdb_sqlite <- args[9]
work_path <- args[10]
get_merged_fasta_R <- args[11]
annotation_R <- args[12]

setwd(work_path)
# Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
message("Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships")
###Load the test data
seurat_with_time <- readRDS(time_rds)
###Get expression profiles ordered by pseudotime
expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 50)
###Filter noise and logFC in expression profile
expression_profile_filter <- filter_expression_profile(expression_profile, FC=0.01)
###K-means clustering
clustering <- clustering_Kmeans(expression_profile_filter, K1=4)
# clustering[1:5,1:5]

Kmeans_clustering_ENS <- clustering
Kmeans_clustering_ENS$Symbol <- rownames(Kmeans_clustering_ENS)
Kmeans_clustering_ENS <- Kmeans_clustering_ENS[, c("Symbol", setdiff(names(Kmeans_clustering_ENS), "Symbol"))]
message("Kmeans clustering result:")
print(Kmeans_clustering_ENS[1:5,1:5])

pdf('plot_kmeans_pheatmap.pdf', width=10, height=10)
plot_kmeans_pheatmap(clustering,ModuleColor1 = c('#67C7C1','#5BA6DA','#FFBF0F','#C067A9'))
dev.off()

message("Loading TFs blinding motif data...")
motif1 <- read.table(tf2motif_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE) # added by yd
message("Checking info of tf blinding motif ")
head(motif1) # added by yd

# # GENIE3[optional] to infer regulatory relationships # Test it for the most small data needed time
# # Ref: https://github.com/aertslab/GENIE3
# library(Seurat)
# library(GENIE3)
# library(doParallel)
# library(doRNG)
# library(reshape2)
# library(dplyr)
# weightMat <- GENIE3(as.matrix(seurat_with_time@assays$RNA@data),nCores = 50) #Error in loadNamespace(x) : there is no package called ‘doParallel’
# weightMat <- getLinkList(weightMat)
# regulation <- weightMat[weightMat[,3]>0.0002,]
# ### add regulation type for each gene pair
# regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS,regulation)
# ### check whether source genes are transcription factors
# motifTF <- c()
# for (i in 1:nrow(motif1)) {
#   TF <- strsplit(motif1[i,5],';')[[1]]
#   motifTF <- c(motifTF,TF)
# }
# regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]
# head(regulatory_relationships)

message("Loading regulatory relationships data from scRNA data runed by GENIE3")
load(genie3_rdata)
head(regulatory_relationships);length(regulatory_relationships$TF)

# # Person’s correlation[optional] to filter regulatory relationships
# regulatory_relationships <- get_cor(Kmeans_clustering_ENS, motif = motif1, correlation_filter = 0.6, start_column = 4)
# head(regulatory_relationships)

# Part 2: Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)
message("Part 2: Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)")
###merge footprints whose distance is less than 4
filtered_footprints <- read.table(footprints_bed,sep = '\t');head(filtered_footprints)
fastadir <- genome_fa

source(get_merged_fasta_R)
merged_fasta <- get_merged_fasta(filtered_footprints,fastadir); head(merged_fasta)
write.table(merged_fasta,'merged_footprints.fasta',row.names=F,quote=F,col.names=F)

### Identify differentially expressed genes related motifs
motif1 <- motifs_select(motif1, rownames(Kmeans_clustering_ENS)) ###Kmeans_clustering_ENS was obtained in part1
### run find_motifs()
fimodir <- 'fimo'
outputdir1 <- './fimo/output/'
outputdir <- './fimo/output/'
motifdir <- pwm_dir
sequencedir <- 'merged_footprints.fasta'
find_motifs(motif1,step=20,fimodir, outputdir1, outputdir, motifdir, sequencedir)

# 文件路径
file_path <- paste0(outputdir1,"Fimo_All.sh")
lines <- readLines(file_path)
last_line <- length(lines)
lines[last_line] <- gsub("&$", "", lines[last_line]) # delete the end &
writeLines(lines, file_path)

### Build fimo env
# # 定义要执行的 shell 命令
# commands <- '
# cd ~
# sudo cp /data/work/SCPipelines/meme-5.5.8.tar.gz .
# tar zxf meme-5.5.8.tar.gz
# cd meme-5.5.8
# ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
# make
# make test
# make install

# export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.8:$PATH
# export PATH=$PATH:~/meme/bin
# '
# # 将命令写入一个临时脚本文件
# temp_script <- tempfile()
# writeLines(commands, temp_script)
# # 执行临时脚本文件
# system(paste("bash", temp_script, ">build_fimo.log 2>&1"))
# # 删除临时脚本文件
# unlink(temp_script)

# 设置环境变量 PATH
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/meme/bin", sep = ":"))
# 执行 fimo --version 并捕获输出
output <- system("fimo --version", intern = TRUE); print(output)

### run fimo_all script in shell
shell_code <- paste0('sh ',outputdir1,'Fimo_All.sh') # delete the end of &
system(shell_code,wait=TRUE)
### delete shell scripts
shell_code2 <- paste0('rm ',outputdir1,'Fimo*.sh')
system(shell_code2,wait=TRUE)

###Combine all footprints of motifs
combined <- combine_footprints(outputdir)
# peaks <- read.delim('differential_peaks.bed')
peaks <- read.delim(peaks_bed)
overlapped <- overlap_footprints_peaks(combined,peaks)
length(overlapped$sequence_name)

###get footprint-related genes
txdb <- loadDb(txdb_sqlite)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # Build species specify TxDb
# sudo /opt/software/miniconda3/envs/IReNA/bin/Rscript -e 'install.packages("/data/work/SCPipelines/build_orgdb/orgdb/org.Ga.eg.db", repos = NULL, type = "sources")'
annodb <- "org.Ga.eg.db"
library(org.Ga.eg.db)

source(annotation_R)
list1 <- get_related_genes(overlapped,txdb = txdb,motif=motif1,Species = 'ga')
str(list1)
###Get candidate genes/TFs-related peaks
list2 <- get_related_peaks(list1,Kmeans_clustering_ENS)
str(list2)
### output filtered footprints
write.table(list2[[1]],'filtered_footprints.bed', quote = F, row.names = F, col.names = F, sep = '\t')

# Define BAM input and output file names
# bam_inputs <- c(
#     "/data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-2_bwa_rmdup.bam",
#     "/data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-1_bwa_rmdup.bam",
#     "/data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1_0_bwa_rmdup.bam",
#     "/data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1_1_bwa_rmdup.bam"
# )
bam_outputs <- paste0("input_", seq_along(bam_inputs), "_filter.bam")

# Run samtools in a loop
for (i in seq_along(bam_inputs)) {
    shell_cmd <- paste(
        "/opt/software/miniconda3/envs/samtools/bin/samtools view -hb -L",
        "filtered_footprints.bed",
        bam_inputs[i],
        ">",
        bam_outputs[i]
    )
    system(shell_cmd, wait = TRUE)
}

# # [optional]
# library(ATACseqQC)
# library(Rsamtools)
# bamfilepath1 <- 'SSC1_filter.bam'
# bamfilepath2 <- 'SSC2_filter.bam'
# bamfilepath3 <- 'esc_filter.bam'
# indexBam(bamfilepath1)
# gal1 <- readBamFile(bamfilepath1, tag=tags, asMates=TRUE, bigFile=TRUE)
# gal2 <- readBamFile(bamfilepath2, tag=tags, asMates=TRUE, bigFile=TRUE)
# gal3 <- readBamFile(bamfilepath3, tag=tags, asMates=TRUE, bigFile=TRUE)
# galout1 <- shiftGAlignmentsList(gal, 'SSC1_filter_shift.bam')
# galout2 <- shiftGAlignmentsList(gal, 'SSC2_filter_shift.bam')
# galout3 <- shiftGAlignmentsList(gal, 'esc_filter_shift.bam')

### calculate cuts of each each position in footprints
bam_filepaths <- bam_outputs
cut_list <- vector("list", length(bam_filepaths))
for (i in seq_along(bam_filepaths)) {
    cut_list[[i]] <- cal_footprint_cuts(
        bamfilepath = bam_filepaths[i],
        bedfile = list2[[1]],
        workers = 8,
        index_bam = TRUE
    )
}

### calculate cuts of each each position in footprints
bamfilepath1 <- "input_1_filter.bam"  
bamfilepath2 <- "input_2_filter.bam"
bamfilepath3 <- "input_3_filter.bam"
bamfilepath4 <- "input_4_filter.bam"
### set parameter 'workers' to make this function run in parallel
cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1,bedfile = list2[[1]],workers = 16,index_bam = T)
cuts2 <- cal_footprint_cuts(bamfilepath = bamfilepath2,bedfile = list2[[1]],workers = 16,index_bam = T)
cuts3 <- cal_footprint_cuts(bamfilepath = bamfilepath3,bedfile = list2[[1]],workers = 16,index_bam = T)
cuts4 <- cal_footprint_cuts(bamfilepath = bamfilepath4,bedfile = list2[[1]],workers = 16,index_bam = T)
cut_list <- list(cuts1,cuts2,cuts3,cuts4)

### get related genes of footprints with high FOS
potential_regulation <- Footprints_FOS(cut_list,list2[[2]], FOS_threshold = 0.1); length(potential_regulation$TF)
# potential_regulation <- Footprints_FOS(cut_list,list2[[2]], FOS_threshold = 0); length(potential_regulation$TF)
### Use information of footprints with high FOS to refine regulatory relationships
# filtered_regulatory <- filter_ATAC(potential_regulation,regulatory_relationships)


filter_ATAC <- function (FOS, regulary_relationships) 
{
    # validInput(regulary_relationships, "regulary_relationships", "df")
    # if (grepl("ENS", FOS[1, 1])) {
    #     TfIndex <- 1
    # } else {
    #     TfIndex <- 2
    # }
    # if (grepl("ENS", FOS[1, 2])) {
    #     TargetIndex <- 4
    # } else {
    #     TargetIndex <- 5
    # }
    TfIndex <- 1; TargetIndex <- 4 # added by yd
    pair1 <- paste(FOS[, 1], FOS[, 2])
    pair2 <- paste(regulary_relationships[, TfIndex], regulary_relationships[, 
        TargetIndex])
    filtered <- regulary_relationships[pair2 %in% pair1, ]
    return(filtered)
}
head(potential_regulation, n=2)
head(regulatory_relationships, n=2)
filtered_regulatory <- filter_ATAC(potential_regulation,regulatory_relationships); head(filtered_regulatory)
length(filtered_regulatory$TF)

TFs_list <- network_analysis(filtered_regulatory,Kmeans_clustering_ENS,TFFDR1 = 10,TFFDR2 = 10, ModuleFDR = 0.05)
str(TFs_list)

save(TFs_list, file = "TFs_list.RData")
# load("TFs_list.RData")
# str(TFs_list)

# Part 3: Regulatory network analysis and visualization
# filtered_regulatory_relationships <- filtered_regulatory # added by yd
# filtered_regulatory_relationships <- regulatory_relationships # added by yd
# TFs_list <- network_analysis(filtered_regulatory_relationships,Kmeans_clustering_ENS,TFFDR1 = 10,TFFDR2 = 10)
# str(TFs_list)

pdf("network.pdf", width=16, height=8)
tryCatch(
    plot_tf_network(TFs_list),
    error = function(e) {
        message("Error in plot_tf_network: ", e$message)
    }
)
dev.off()

# load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
# Kmeans_cluster_Ens <- add_ENSID(test_clustering,Spec1='Hs')
# head(Kmeans_cluster_Ens)
# motif1 <- Tranfac201803_Hs_MotifTFsF
# head(motif1)
# regulatory_relationships <- get_cor(Kmeans_cluster_Ens,motif1,0.7); head(regulatory_relationships)
# TFs_list <- network_analysis(regulatory_relationships,Kmeans_cluster_Ens)
# str(TFs_list)

### Download Homo sapiens org.db
#BiocManger::install('org.Hs.eg.db')
# library(org.Ga.eg.db)
# ### Enrichment analysis (KEGG)
# enrichment_KEGG <- enrich_module(Kmeans_clustering_ENS, org.Ga.eg.db, enrich.db = 'KEGG',organism = 'ga')
# #enrichment_GO <- enrich_module(Kmeans_cluster_ENS, org.Hs.eg.db, 'GO')

# ### Enrichment analysis (GO)
# enrichment_GO <- enrich_module(Kmeans_clustering_ENS, enrich.db = 'GO',org.Ga.eg.db)