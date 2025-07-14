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

# # Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
genie3_rdata="/data/work/test0711/IReNA_GENIE3_2/seu_day-2.rds_seurat_with_time.rds_GENIE3.RData"

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

# Part 2: Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)
filtered_footprints_bed='/data/work/test_raw/10_footprints/filtered_footprints.bed'
fastadir="~/rgtdata/tair10/genome_tair10_ensembl_release_51.fa"
get_merged_fasta_R="/data/work/SCPipelines/all/get_merged_fasta.R"
motifdir='/data/work/SCPipelines/Gar_individual_motif2/'
peaks_bed='/data/work/test_raw/08_differential_peaks/peaks.bed'
txdb_sqlite="/data/work/SCPipelines/gtf2txdb/txdb.sqlite"
annotation_R="/data/work/SCPipelines/all/annotation.R"
n_workers=8
annodb <- "org.Ga.eg.db"

###merge footprints whose distance is less than 4
filtered_footprints <- read.table(filtered_footprints_bed,sep = '\t');head(filtered_footprints)
source(get_merged_fasta_R)
merged_fasta <- get_merged_fasta(filtered_footprints,fastadir); head(merged_fasta)
write.table(merged_fasta,'merged_footprints.fasta',row.names=F,quote=F,col.names=F)

### Identify differentially expressed genes related motifs
find_and_overlap_footprints <- function(
    motif1,
    kmeans_clustering_ens,
    motifdir,
    peaks_bed,
    sequencedir = "merged_footprints.fasta",
    fimodir = "fimo",
    outputdir1 = "./fimo/output/",
    outputdir = "./fimo/output/"
) {
    motif1 <- motifs_select(motif1, rownames(kmeans_clustering_ens))

    if (!dir.exists(fimodir)) {
        dir.create(fimodir, recursive = TRUE)
    }
    if (!dir.exists(outputdir1)) {
        dir.create(outputdir1, recursive = TRUE)
    }

    find_motifs(motif1,step=20,fimodir, outputdir1, outputdir, motifdir, sequencedir)

    file_path <- paste0(outputdir1, "Fimo_All.sh")
    lines <- readLines(file_path)
    last_line <- length(lines)
    lines[last_line] <- gsub("&$", "", lines[last_line])
    writeLines(lines, file_path)

    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/meme/bin", sep = ":"))
    output <- system("fimo --version", intern = TRUE)
    print(output)

    shell_code <- paste0("sh ", outputdir1, "Fimo_All.sh")
    system(shell_code, wait = TRUE)

    shell_code2 <- paste0("rm ", outputdir1, "Fimo*.sh")
    system(shell_code2, wait = TRUE)

    combined <- combine_footprints(outputdir)
    peaks <- read.delim(peaks_bed)
    save(combined, peaks, file = "combined_peaks.RData")
    # overlapped <- overlap_footprints_peaks(combined, peaks)
    # return(overlapped)
}
find_and_overlap_footprints(
    motif1 = motif1,
    kmeans_clustering_ens = Kmeans_clustering_ENS,
    motifdir = motifdir,
    peaks_bed = peaks_bed,
)


# overlap_footprints_peaks <- function (footprints, peak_bed) 
# {
#     # validInput(footprints, "footprints", "df")
#     # validInput(peak_bed, "peak_bed", "df")
#     peak <- GenomicRanges::GRanges(paste0(peak_bed[, 1], ":", 
#         peak_bed[, 2], "-", peak_bed[, 3]))
#     message("First step")
#     footprints1 <- GenomicRanges::GRanges(paste0(footprints[, 
#         1], ":", footprints[, 2], "-", footprints[, 3]))
#     message("Second step")
#     overlap_idx <- GenomicRanges::findOverlaps(footprints1, peak)
#     message("Third step")
#     final_footprints <- cbind(footprints[overlap_idx@from, ], 
#         peak_bed[overlap_idx@to, ])
#     message("End step")
#     return(final_footprints)
# }
# overlapped <- overlap_footprints_peaks(combined,peaks)
############################## overlap_footprints_peaks ############################
# footprints <- combined
# peak_bed <- peaks
# peak <- GenomicRanges::GRanges(paste0(peak_bed[, 1], ":", peak_bed[, 2], "-", peak_bed[, 3]))
# message("First step")
# footprints1 <- GenomicRanges::GRanges(paste0(footprints[, 1], ":", footprints[, 2], "-", footprints[, 3]))
# message("Second step")
# overlap_idx <- GenomicRanges::findOverlaps(footprints1, peak)
# message("Third step")
# final_footprints <- cbind(footprints[overlap_idx@from, ], peak_bed[overlap_idx@to, ])
# message("End step")
# overlapped <- final_footprints
# length(overlapped$sequence_name)
# message("passed overlapping")

###get footprint-related genes
txdb <- loadDb(txdb_sqlite)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # Build species specify TxDb
# sudo /opt/software/miniconda3/envs/IReNA/bin/Rscript -e 'install.packages("/data/work/SCPipelines/build_orgdb/orgdb/org.Ga.eg.db", repos = NULL, type = "sources")'
annodb <- "org.Ga.eg.db"
do.call(library, list(annodb))
source(annotation_R)
list1 <- get_related_genes(overlapped,txdb = txdb,motif=motif1,Species = 'ga')
str(list1)
###Get candidate genes/TFs-related peaks
list2 <- get_related_peaks(list1,Kmeans_clustering_ENS)
str(list2)
### output filtered footprints
write.table(list2[[1]],'filtered_footprints.bed', quote = F, row.names = F, col.names = F, sep = '\t')

# ### run samtools in shell
bam_inputs <- c('/data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-2_bwa_rmdup.bam')
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
bam_filepaths <- bam_outputs
cut_list <- vector("list", length(bam_filepaths))
for (i in seq_along(bam_filepaths)) {
    cut_list[[i]] <- cal_footprint_cuts(
        bamfilepath = bam_filepaths[i],
        bedfile = list2[[1]],
        workers = n_workers,
        index_bam = T
    )
}
# bamfilepath1 <- "input_1_filter.bam"  
# ### set parameter 'workers' to make this function run in parallel
# cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1,bedfile = list2[[1]],workers = n_workers,index_bam = T)
# cut_list <- list(cuts1)