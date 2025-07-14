# IReNA /opt/software/miniconda3/envs/IReNA/bin/R
library(IReNA)
library(AnnotationDbi)
library(Seurat)
library(GENIE3)
library(doParallel)
library(doRNG)
library(reshape2)
library(dplyr)
library(clusterProfiler)
library(optparse)

option_list <- list(
    make_option("--genie3_rdata", type = "character", default = "/data/work/N2/ppscRNA/seu_day-2.rds_hvg3000.rds_GENIE3.RData", help = "Path to GENIE3 RData file"),
    make_option("--filtered_footprints_bed", type = "character", default = "/data/work/N2/10_footprints/filtered_footprints.bed", help = "Filtered footprints BED file"),
    make_option("--fastadir", type = "character", default = "~/rgtdata/tair10/genome_tair10_ensembl_release_51.fa", help = "FASTA directory"),
    make_option("--get_merged_fasta_R", type = "character", default = "/data/work/SCPipelines/all/get_merged_fasta.R", help = "get_merged_fasta R script"),
    make_option("--motifdir", type = "character", default = "/data/work/SCPipelines/Gar_individual_motif2/", help = "Motif directory"),
    make_option("--peaks_bed", type = "character", default = "/data/work/N2/08_differential_peaks/peaks.bed", help = "Peaks BED file"),
    make_option("--txdb_sqlite", type = "character", default = "/data/work/SCPipelines/gtf2txdb/txdb.sqlite", help = "TxDb sqlite file"),
    make_option("--annotation_R", type = "character", default = "/data/work/SCPipelines/all/annotation.R", help = "Annotation R script"),
    make_option("--n_workers", type = "integer", default = 8, help = "Number of workers"),
    make_option("--tf_fdr", type = "integer", default = 2, help = "TF FDR"),
    make_option("--annodb", type = "character", default = "org.Ga.eg.db", help = "Annotation DB"),
    make_option("--bams_txt", type = "character", default = "/data/work/N2/bams.txt", help = "BAMs txt file"),
    make_option("--work_path", type = "character", default = "/data/work/N2/bulkATAC_scRNA", help = "Working path")
)

opt <- parse_args(OptionParser(option_list = option_list))

genie3_rdata <- opt$genie3_rdata
filtered_footprints_bed <- opt$filtered_footprints_bed
fastadir <- opt$fastadir
get_merged_fasta_R <- opt$get_merged_fasta_R
motifdir <- opt$motifdir
peaks_bed <- opt$peaks_bed
txdb_sqlite <- opt$txdb_sqlite
annotation_R <- opt$annotation_R
n_workers <- opt$n_workers
tf_fdr <- opt$tf_fdr
annodb <- opt$annodb
bams_txt <- opt$bams_txt
work_path <- opt$work_path

setwd(work_path)

# # Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
# ###Load the test data
# seurat_with_time <- readRDS('/data/work/test0711/seu_day-2.rds_seurat_with_time.rds')
# ###Get expression profiles ordered by pseudotime
# expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 50)
# ###Filter noise and logFC in expression profile
# expression_profile_filter <- filter_expression_profile(expression_profile, FC=0.01)
# ###K-means clustering
# clustering <- clustering_Kmeans(expression_profile_filter, K1=4)
# # clustering[1:5,1:5]

# Kmeans_clustering_ENS <- clustering
# Kmeans_clustering_ENS$Symbol <- rownames(Kmeans_clustering_ENS)
# Kmeans_clustering_ENS <- Kmeans_clustering_ENS[, c("Symbol", setdiff(names(Kmeans_clustering_ENS), "Symbol"))]
# Kmeans_clustering_ENS[1:5,1:5]
# length(unique(rownames(Kmeans_clustering_ENS)))

# pdf('plot_kmeans_pheatmap.pdf', width=10, height=10)
# plot_kmeans_pheatmap(clustering,ModuleColor1 = c('#67C7C1','#5BA6DA','#FFBF0F','#C067A9'))
# dev.off()

# motif_txt <- "/data/work/SCPipelines/pp_scRNA/tf_blinding_motif_ga.txt"
# motif1 <- read.table(motif_txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE) # added by yd
# head(motif1) # added by yd
# length(unique(motif1$EnsemblID))

# # # GENIE3[optional] to infer regulatory relationships # Test it for the most small data needed time
# # # Ref: https://github.com/aertslab/GENIE3
# # weightMat <- GENIE3(as.matrix(seurat_with_time@assays$RNA@data),nCores = 50) #Error in loadNamespace(x) : there is no package called ‘doParallel’
# # weightMat <- getLinkList(weightMat)
# # regulation <- weightMat[weightMat[,3]>0.0002,]
# # ### add regulation type for each gene pair
# # regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS,regulation)
# # ### check whether source genes are transcription factors
# # motifTF <- c()
# # for (i in 1:nrow(motif1)) {
# #   TF <- strsplit(motif1[i,5],';')[[1]]
# #   motifTF <- c(motifTF,TF)
# # }
# # regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]
# # head(regulatory_relationships)

load(genie3_rdata)
message("load regulatory_relationships and print 2 rows:")
print(head(regulatory_relationships,n=2))
length(regulatory_relationships$TF)
message("load Kmeans_clustering_ENS and print 2 rows:")
print(head(Kmeans_clustering_ENS,n=2))
length(Kmeans_clustering_ENS$Symbol)
message("load motif1 and print 2 rows:")
print(head(motif1,n=2))
length(motif1$TFs)

# # Person’s correlation[optional] to filter regulatory relationships
# regulatory_relationships <- get_cor(Kmeans_clustering_ENS, motif = motif1, correlation_filter = 0.6, start_column = 4)
# head(regulatory_relationships)

# Part 2: Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)
###merge footprints whose distance is less than 4
filtered_footprints <- read.table(filtered_footprints_bed,sep = '\t');head(filtered_footprints)
source(get_merged_fasta_R)
merged_fasta <- get_merged_fasta(filtered_footprints,fastadir); head(merged_fasta)
write.table(merged_fasta,'merged_footprints.fasta',row.names=F,quote=F,col.names=F)

### Identify differentially expressed genes related motifs
motif1 <- motifs_select(motif1, rownames(Kmeans_clustering_ENS)) ###Kmeans_clustering_ENS was obtained in part1 
length(motif1$TFs)
### run find_motifs()
fimodir <- 'fimo'
outputdir1 <- './fimo/output/'
outputdir <- './fimo/output/'
motifdir <- motifdir
sequencedir <- 'merged_footprints.fasta'
if (!dir.exists(fimodir)) {
  dir.create(fimodir, recursive = TRUE)
}
if (!dir.exists(outputdir1)) {
  dir.create(outputdir1, recursive = TRUE)
}
find_motifs(motif1,step=20,fimodir, outputdir1, outputdir, motifdir, sequencedir)

# 文件路径
# file_path <- paste0(outputdir1,"Fimo_All.sh")
# lines <- readLines(file_path)
# last_line <- length(lines)
# lines[last_line] <- gsub("&$", "", lines[last_line]) # delete the end &
# writeLines(lines, file_path)

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
peaks <- read.delim(peaks_bed,header = FALSE)
# peaks <- read.delim(peaks_bed)
# overlapped <- overlap_footprints_peaks(combined,peaks)

# 定义最大尝试次数
max_attempts <- 10
attempts <- 0
# 循环尝试运行代码，直到成功或达到最大尝试次数
while (attempts < max_attempts) {
  attempts <- attempts + 1
  tryCatch({
    overlapped <- overlap_footprints_peaks(combined, peaks)
    message("Success! The code ran without errors.")
    break  # 如果成功，退出循环
  }, error = function(e) {
    message(paste("Attempt", attempts, "failed with error:", e$message))
    Sys.sleep(60)  # 等待一段时间后再次尝试，避免过快重复尝试
  })
}

# 检查是否成功
if (attempts == max_attempts) {
  message("The code failed after", max_attempts, "attempts. Please check for issues.")
} else {
  message("The code ran successfully after", attempts, "attempts.")
}

###get footprint-related genes
txdb <- loadDb(txdb_sqlite)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # Build species specify TxDb
# sudo /opt/software/miniconda3/envs/IReNA/bin/Rscript -e 'install.packages("/data/work/SCPipelines/build_orgdb/orgdb/org.Ga.eg.db", repos = NULL, type = "sources")'
do.call(library, list(annodb))
source(annotation_R)

# 定义最大尝试次数
max_attempts <- 10
attempts <- 0
# 循环尝试运行代码，直到成功或达到最大尝试次数
while (attempts < max_attempts) {
  attempts <- attempts + 1
  tryCatch({
    list1 <- get_related_genes(overlapped,txdb = txdb,motif=motif1,annodb = annodb)
    message("Success! The code ran without errors.")
    break  # 如果成功，退出循环
  }, error = function(e) {
    message(paste("Attempt", attempts, "failed with error:", e$message))
    Sys.sleep(60)  # 等待一段时间后再次尝试，避免过快重复尝试
  })
}

# 检查是否成功
if (attempts == max_attempts) {
  message("The code failed after", max_attempts, "attempts. Please check for issues.")
} else {
  message("The code ran successfully after", attempts, "attempts.")
}


str(list1)
###Get candidate genes/TFs-related peaks
list2 <- get_related_peaks(list1,Kmeans_clustering_ENS)
str(list2)
### output filtered footprints
write.table(list2[[1]],'filtered_footprints.bed', quote = F, row.names = F, col.names = F, sep = '\t')

# ### run samtools in shell
# 读取文件内容
bam_inputs <- readLines(bams_txt)
bam_inputs <- strsplit(bam_inputs, ",")[[1]]
print(bam_inputs)

bam_outputs <- paste0("input_", seq_along(bam_inputs), "_filter.bam")
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

# shell_code1 <- '/opt/software/miniconda3/envs/samtools/bin/samtools view -hb -L filtered_footprints.bed /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-2_bwa_rmdup.bam > E1_N2_filter.bam'
# shell_code2 <- '/opt/software/miniconda3/envs/samtools/bin/samtools view -hb -L filtered_footprints.bed /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-1_bwa_rmdup.bam > E1_N1_filter.bam'
# shell_code3 <- '/opt/software/miniconda3/envs/samtools/bin/samtools view -hb -L filtered_footprints.bed /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1_0_bwa_rmdup.bam > E1_0_filter.bam'
# shell_code4 <- '/opt/software/miniconda3/envs/samtools/bin/samtools view -hb -L filtered_footprints.bed /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1_1_bwa_rmdup.bam > E1_P1_filter.bam'
# bam_outputs <- paste0("input_", seq_along(bam_inputs), "_filter.bam")
# system(shell_code1,,wait=TRUE)
# system(shell_code2,,wait=TRUE)
# system(shell_code3,,wait=TRUE)
# system(shell_code4,,wait=TRUE)

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

# ### calculate cuts of each each position in footprints
# 定义最大尝试次数
max_attempts <- 10
attempts <- 0
# 循环尝试运行代码，直到成功或达到最大尝试次数
while (attempts < max_attempts) {
  attempts <- attempts + 1
  tryCatch({
    cut_list <- vector("list", length(bam_outputs))
    for (i in seq_along(bam_outputs)) {
        cut_list[[i]] <- cal_footprint_cuts(
            bamfilepath = bam_outputs[i],
            bedfile = list2[[1]],
            workers = n_workers,
            index_bam = T
        )
    }
    message("Success! The code ran without errors.")
    break  # 如果成功，退出循环
  }, error = function(e) {
    message(paste("Attempt", attempts, "failed with error:", e$message))
    Sys.sleep(60)  # 等待一段时间后再次尝试，避免过快重复尝试
  })
}

# 检查是否成功
if (attempts == max_attempts) {
  message("The code failed after", max_attempts, "attempts. Please check for issues.")
} else {
  message("The code ran successfully after", attempts, "attempts.")
}


# bamfilepath1 <- "input_1_filter.bam"  
# ### set parameter 'workers' to make this function run in parallel
# cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1,bedfile = list2[[1]],workers = n_workers,index_bam = T)
# cut_list <- list(cuts1)
# bamfilepath1 <- 'E1_N2_filter.bam'
# bamfilepath2 <- 'E1_N1_filter.bam'
# bamfilepath3 <- 'E1_0_filter.bam'
# bamfilepath4 <- 'E1_P1_filter.bam'
# ### set parameter 'workers' to make this function run in parallel
# cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1,bedfile = list2[[1]],workers = 32,index_bam = T)
# cuts2 <- cal_footprint_cuts(bamfilepath = bamfilepath2,bedfile = list2[[1]],workers = 32,index_bam = T)
# cuts3 <- cal_footprint_cuts(bamfilepath = bamfilepath3,bedfile = list2[[1]],workers = 32,index_bam = T)
# cuts4 <- cal_footprint_cuts(bamfilepath = bamfilepath4,bedfile = list2[[1]],workers = 32,index_bam = T)
# cut_list <- list(cuts1,cuts2,cuts3,cuts4)

### get related genes of footprints with high FOS
potential_regulation <- Footprints_FOS(cut_list,list2[[2]], FOS_threshold = 0.1); length(potential_regulation$TF)
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

TFs_list <- network_analysis(filtered_regulatory,Kmeans_clustering_ENS,TFFDR1 = tf_fdr,TFFDR2 = tf_fdr)
# TFs_list <- network_analysis(filtered_regulatory,Kmeans_clustering_ENS,TFFDR1 = 1,TFFDR2 = 1)
# str(TFs_list)

# Part 3: Regulatory network analysis and visualization
# filtered_regulatory_relationships <- filtered_regulatory # added by yd

pdf("network.pdf", height=8, width=12)
tryCatch({
  plot_tf_network(TFs_list, layout = "grid")
}, error = function(e) {
  message("Error in plotting: ", e$message)
})
dev.off()

save(TFs_list, file = "TFs_list_.RData")
save(filtered_regulatory, file = "filtered_regulatory.RData")
save(potential_regulation, file = "potential_regulation.RData")
save(Kmeans_clustering_ENS,file = "Kmeans_clustering_ENS.RData")
save(motif1,file = "motif1.RData")
message(paste0("Saved results TFs_list, filtered_regulatory, potential_regulation, Kmeans_clustering_ENS, motif1 in: ", basename(genie3_rdata), "_regulatory.RData"))

enrich_module <- function(Kmeans_result, org.db, enrich.db ,fun_num = 5,
                          pvalueCutoff = 0.05, use_internal_data = TRUE, organism = NULL) {
  all_gene <- Kmeans_result
  all_gene<-all_gene[order(all_gene$KmeansGroup),]
  le<-levels(as.factor(all_gene$KmeansGroup))
  for (i in le) {
    acc11<-rownames(all_gene[all_gene$KmeansGroup == i,])
    # gene1 <- clusterProfiler::bitr(acc11, fromType = "ENSEMBL",
    #               toType = c("SYMBOL", "ENTREZID"),
    #               OrgDb = org.db)
    gene1 <- data.frame(ENTREZID = acc11, row.names = acc11)
    if (enrich.db =='KEGG') {
      k1 <- clusterProfiler::enrichKEGG(gene = gene1$ENTREZID,
                                        pvalueCutoff = pvalueCutoff
                                        ,organism = organism
                                        ,use_internal_data = use_internal_data)
      # kegg_result <- enricher(gene_list,TERM2GENE = pathway2gene, TERM2NAME = pathway2name,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    }else if(enrich.db =='GO'){
      k1 = clusterProfiler::enrichGO(gene = gene1$ENTREZID,
                    OrgDb = org.db,
                    keyType = "GID",
                    ont = "BP",
                    pvalueCutoff = pvalueCutoff)
    }
    acc2 <- k1@result
    acc2$'-log10(q-value)' <- -log10(acc2$qvalue)
    acc2 <- acc2[order(-acc2$`-log10(q-value)`),]
    if (i=='1' | i == 1) {
      acc21 <- acc2[1:fun_num,]
      acc21$module <- rep(i,fun_num)
    }else{
      acc22 <- acc2[1:fun_num,]
      acc22$module<-rep(i,fun_num)
      acc21 <- rbind(acc21,acc22)
    }
  }
  acc21 <- acc21[,c(1,2,11,10,3:9)]
  return(acc21)
}
enrichment_GO <- enrich_module(Kmeans_clustering_ENS, org.db=annodb, enrich.db = 'GO', fun_num = 5, pvalueCutoff = 0.05, use_internal_data = TRUE)
### select functions that you want to present in the figure
enrichment_GO_select <- enrichment_GO[c(1,11,21),]
enrichment_GO_select <- na.omit(enrichment_GO_select)
### plotting
pdf('go.pdf', height = 8, width = 12)

tryCatch({
  plot_intramodular_network(TFs_list, enrichment_GO_select, layout = 'circle')
}, error = function(e) {
  message("Error in plotting: ", e$message)
})

dev.off()