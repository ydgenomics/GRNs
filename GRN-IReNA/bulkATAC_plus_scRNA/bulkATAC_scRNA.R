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
##################### get_merged_fasta #####################
get_merged_fasta <- function(footprints, fastadir, distance = 4) {
  # validInput(footprints,'footprints','df') # added by yd
  # validInput(fastadir,'fastadir','direxists') # added by yd
  # validInput(distance,'distance','numeric') # added by yd
  merged_footprints <- merge_footprints(footprints, ditance = distance)
  fasta <- getfasta(merged_footprints, fastadir = fastadir)
  return(fasta)
}


merge_footprints <- function(footprints, ditance = 4, revise = TRUE) {
  con1 <- footprints
  str1 <- con1[1,1]
  start1 <- con1[1,2]
  end1 <- con1[1,3]
  pva1 <- con1[1,4]
  no1 <- 1
  out1 <- c()
  for (i in 2:nrow(con1)) {
    str2 <- as.character(con1[i, ][1])
    start2 <- as.numeric(con1[i, ][2])
    end2 <- as.numeric(con1[i, ][3])
    pva2 <- as.numeric(con1[i, ][4])
    if (str2 == str1 & (start1 >= start2 & start1 <= end2 | end1 >= start2 &
                        end1 <= end2 | start1 <= start2 & end1 >= end2 | start2 >
                        end1 & (start2 - end1) <= ditance | start1 > end2 &
                        (start1 - end2) <= ditance)) {
      acc22 <- sort(c(start1, start2, end1, end2))
      start1 <- acc22[1]
      end1 <- acc22[4]
      pva1 <- pva1 + pva2
      no1 <- no1 + 1
    } else {
      pva2 <- pva1 / no1
      col1 <- paste(str1, start1, end1, pva2, sep = "\t")
      out1 <- c(out1, col1)
      str1 <- str2
      start1 <- start2
      end1 <- end2
      pva1 <- as.numeric(con1[i, ][4])
      no1 <- 1
    }
  }
  out2 <- as.data.frame(t(as.data.frame(strsplit(out1, "\t"))))
  colnames(out2) <- c("chr", "strat", "stop", "pvalue")
  out2[, 2] <- as.numeric(out2[, 2])
  out2[, 3] <- as.numeric(out2[, 3])
  if (revise == TRUE) {
    out2[, 2] <- out2[, 2] - 2
    out2[, 3] <- out2[, 3] + 2
    return(out2)
  } else {
    return(out2)
  }
}


getfasta <- function(merged_footprints, fastadir) {
  fasta <- Biostrings::readBStringSet(fastadir, format = "fasta", nrec = -1L,
                                      skip = 0L, seek.first.rec = FALSE,
                                      use.names = TRUE)
  fasta1 <- c()
  for (i in 1:nrow(merged_footprints)) {
    name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                   "-", merged_footprints[i, 3])
    if (merged_footprints[i, 3] > length(fasta[[merged_footprints[i, 1]]])) {
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):length(fasta[[merged_footprints[i, 1]]])]))
      name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                     "-", length(fasta[[merged_footprints[i, 1]]]))
    } else{
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):merged_footprints[i, 3]]))
    }
    fasta1 <- c(fasta1, name, sequence)
  }
  fasta1 <- as.data.frame(fasta1)
  return(fasta1)
}
################################################################
# source("/data/work/SCPipelines/all/get_merged_fasta.R")
merged_fasta <- get_merged_fasta(filtered_footprints,fastadir); head(merged_fasta)
write.table(merged_fasta,'merged_footprints.fasta',row.names=F,quote=F,col.names=F)

### Identify differentially expressed genes related motifs
motif1 <- motifs_select(motif1, rownames(Kmeans_clustering_ENS)) ###Kmeans_clustering_ENS was obtained in part1
### run find_motifs()
fimodir <- 'fimo'
outputdir1 <- './fimo/output/'
outputdir <- './fimo/output/'
motifdir <- '/data/work/SCPipelines/Gar_individual_motif2/'
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
txdb <- loadDb("/data/work/SCPipelines/gtf2txdb/txdb.sqlite")
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # Build species specify TxDb
# sudo /opt/software/miniconda3/envs/IReNA/bin/Rscript -e 'install.packages("/data/work/SCPipelines/build_orgdb/orgdb/org.Ga.eg.db", repos = NULL, type = "sources")'
annodb <- "org.Ga.eg.db"
library(org.Ga.eg.db)
# source("/data/work/SCPipelines/all/annotation.R")
############################### annotation.R ###############################
# Ref: https://rdrr.io/github/jiang-junyao/IReNA/src/R/annotation.R annotation.R
#' Annotate peaks based on regions of peak
#' @description This function first merge and extend footprint regions, and then
#'  integrate R package ChIPseeker to get footprint-related genes
#' @param footprints footprints that overlap with peaks, generated by overlap_footprints_peaks() or
#'  intersect function of bedtools
#' @param motif motif file, you can choose our bulit-in motif database of
#' 'mus musculus', 'homo sapiens', 'zebrafish' and 'chicken' by 'motif = Tranfac201803_Mm_MotifTFsF',
#' 'motif = Tranfac201803_Hs_MotifTFsF', 'motif = Tranfac201803_Zf_MotifTFsF',
#' 'motif = Tranfac201803_Ch_MotifTFsF' respectively, or you can upload your own motif data base,
#'  but the formata use be the same as our built-in motif database.
#' @param Species character, indicating the species of data which is used to
#' choose annodb when annotate peak.
#' @param txdb 	TxDb object contained transcript-related features of a particular
#' genome. Bioconductor provides several package that containing TxDb object of
#' model organisms with multiple commonly used genome version, for instance
#' TxDb.Hsapiens.UCSC.hg38.knownGene, TxDb.Hsapiens.UCSC.hg19.knownGene for
#' human genome hg38 and hg19, TxDb.Mmusculus.UCSC.mm10.knownGene and TxDb.Mmusculus.UCSC.mm9.knownGene
#' for mouse genome mm10 and mm9, etc.
#' @param tssRegion Region Range of TSS
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom ChIPseeker annotatePeak
#' @return return a list, first element is bed format datafrmae, second element
#' is annotated footprints dataframe
#' @export
#'
#' @examples #txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
#' load(system.file("extdata", "combined.rda", package = "IReNA"))
#' load(system.file("extdata", "test_peak.rda", package = "IReNA"))
#' peak_bed <- get_bed(test_peak)
#' overlapped <- overlap_footprints_peaks(combined, peak_bed)
#' #list1 <- get_related_genes(overlapped,txdb = txdb,motif=Tranfac201803_Mm_MotifTFsF,Species = 'Mm')

# footprints <- overlapped
# motif <- motif1
# txdb <- txdb
# Species <- "ga"
# tssRegion = c(-3000, 3000)

get_related_genes <- function(footprints, motif, Species, txdb, tssRegion = c(-3000, 3000)) {
  # validInput(footprints,'footprints','df')
  # validInput(motif,'motif','df')
  # validInput(Species,'Species','character')
  # validInput(txdb,'txdb','txdb')
  # validInput(tssRegion,'tssRegion','vector')
  footprintslist <- merge_extent_footprints(footprints, motif)
  merged_footprints <- footprintslist[[2]]
  if (Species == "Hs") {
    annodb <- "org.Hs.eg.db"
  } else if (Species == "Mm") {
    annodb <- "org.Mm.eg.db"
  } else if (Species == "Zf") {
    annodb <- "org.Dr.eg.db"
  } else if (Species == "Ch") {
    annodb <- "org.Gg.eg.db"
  } else if (Species == "ga") { # added by yd
    annodb <- "org.Ga.eg.db"
  }
  reference_GRange <- GenomicRanges::GRanges(seqnames = merged_footprints[,1],
                                             IRanges::IRanges(start = as.numeric(merged_footprints[,2]),
                                                              end = as.numeric(merged_footprints[,3])),
                                             strand = merged_footprints[,4])
  peakAnno <- ChIPseeker::annotatePeak(reference_GRange,
    tssRegion = tssRegion,
    TxDb = txdb, annoDb = annodb
  )
  region <- peakAnno@anno@elementMetadata$annotation
  # gene <- peakAnno@anno@elementMetadata$ENSEMBL
  gene <- peakAnno@anno@elementMetadata$geneId # added by yd
  start1 <- peakAnno@anno@ranges@start
  merged_footprints2 <- merged_footprints[merged_footprints$V2 %in% start1, ]
  exon1 <- grep('exon',region)
  Intron1 <- grep('Intron',region)
  Intergenic1 <- grep('Intergenic',region)
  Downstream1 <- grep('Downstream',region)
  Promoter1 <- grep('Promoter',region)
  UTR3 <- grep("3' UTR",region)
  UTR5 <- grep("5' UTR",region)
  region2 <- rep(NA,length(region))
  region2[exon1]='Exon'
  region2[Intron1]='Intron'
  region2[Downstream1]='Downstream'
  region2[Promoter1]='Promoter'
  region2[UTR3]="3' UTR"
  region2[UTR5]="5' UTR"
  region2[Intergenic1]='Intergenic'
  table(region2)
  peak_region1 <- paste(as.character(peakAnno@anno@seqnames),
                        as.character(peakAnno@anno@ranges),sep = ':')
  peak_region2 <- paste0(merged_footprints[,1],':'
                         ,merged_footprints[,2],'-',merged_footprints[,3])
  merged_footprints2 <- merged_footprints[peak_region2%in%peak_region1,]
  merged_footprints2$gene <- gene
  merged_footprints2$region <- region2
  merged_footprints2 <- merged_footprints2[, c(9, 8, 1:7)]
  colnames(merged_footprints2) <- c(paste0("V", 1:9))
  footprintslist[[2]] <- merged_footprints2
  footprintslist[[1]] <- footprintslist[[1]][peak_region2%in%peak_region1,]
  return(footprintslist)
}



merge_extent_footprints <- function(file1, motif1) {
  file1 <- file1[file1[,7] %in% motif1$Accession,]
  peak_region <- paste(as.character(file1[,1]), file1[,2],
                       file1[,3], file1[,4], sep = "\t")
  file1$peak_region <- peak_region
  file1$tf <- motif1[match(file1[,7],motif1[,1]),5]
  group_peak <- dplyr::group_by(file1,peak_region)
  group_peak <- group_peak[order(group_peak$peak_region),]
  group_peak2 <-dplyr::group_map(group_peak,~merge_group(.x))
  group_peak3 <- as.data.frame(group_peak[!duplicated(group_peak$peak_region),])
  num1 <- as.integer((as.numeric(group_peak3[,2]) + as.numeric(group_peak3[,3])) / 2)
  size1 <- as.numeric(group_peak3[,3]) - as.numeric(group_peak3[,2]) + 1
  num11 <- as.integer(num1 - (size1 * 5))
  num12 <- as.integer(num1 + (size1 * 5))
  peak_index = paste0(group_peak3[,8],':',group_peak3[,9],'-',group_peak3[,10])
  Candid <- group_peak3[,1:4]
  Candid <- dplyr::mutate(Candid, V5 = num11 ,V6=num12,V7=unlist(group_peak2))
  bed <- Candid[,c(1,5,6)]
  list1 <- list(bed,Candid)
  return(list1)
}



merge_group <- function(data1){
  data1 <- as.data.frame(data1)
  motif <- as.character(data1[,7])
  related_gene <- as.character(data1[,11])
  mg <- paste(motif,related_gene,sep = ';',collapse = '|')
  return(mg)
}
##########################################################
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