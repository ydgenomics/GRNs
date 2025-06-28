# Title: extract_promoters.R
# Date: 20250526
# Coder: ydgenomics
# Description:
# Input: gtf file and fasta file
# Output: 3kpromoter.fasta
# Image: GRN-SCENIC-database--01 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ; https://github.com/aertslab/create_cisTarget_databases

library(rtracklayer)
library(tidyverse)
library(Biostrings)
library(optparse)

option_list <- list(
  make_option(c("-g", "--gtf"),type = "character",default = "/data/work/0.peanut/GRN/AT/Athaliana_TAIR10.54.gtf",help = "Input gtf file"),
  make_option(c("-f", "--fasta"),type = "character",default = "/data/work/0.peanut/GRN/AT/Athaliana_TAIR10.54.dna.fa",help = "Input FASTA file"),
  make_option(c("-s", "--species"),type = "character",default = "AT",help = "Species name, e.g., peanut")
)

opt <- parse_args(OptionParser(option_list = option_list))
input_gtf <- opt$gtf
input_fasta <- opt$fasta
species <- opt$species

gtf <- import(input_gtf, format = "gtf") %>%
  as.data.frame()

head(gtf)
colnames(gtf)
# [1] "seqnames"      "start"         "end"           "width"        
# [5] "strand"        "source"        "type"          "score"        
# [9] "phase"         "transcript_id" "gene_id"
unique(gtf$type)
#[1] transcript exon       CDS

# get gene promoters regions
genes <- gtf %>%
  dplyr::filter(type == "CDS" & gene_id != "NA") %>%
  dplyr::select(seqnames, start, end, strand, gene_id) %>%
  dplyr::mutate(start2 = ifelse(strand == "+", start - 3000, end + 1),
                end2 = ifelse(strand == "+", start - 1, end + 3000))
#对于正链（strand == "+"），启动子区域从 start - 3000 到 start - 1。
#对于负链（strand == "-"），启动子区域从 end + 1 到 end + 3000。
genes$seqnames <- as.character(genes$seqnames)
#genes$seqnames <- sapply(strsplit(genes$seqnames, split = "\\."), "[", 4)
#head(genes$seqnames, 5)
head(genes) # genes$seqnames必须要和names(genome)一致

########## Rename chromosomes of genome ##########
# load genome
genome <- readDNAStringSet(input_fasta)
genome

# rename seqnames
names(genome) <- sapply(strsplit(names(genome),split = " "), "[", 1) #unmatch peanut genome
#names(genome) <- sapply(strsplit(names(genome), split = "\\."), "[", 4)
head(names(genome), 25)
unique(names(genome))
unique(genes$seqnames)

########### Cycle extracting sequences of promoter ###########
# extract promter sequences
# x = 1
lapply(1:nrow(genes), function(x) {
  print(x)
  tmp <- genes[x, ]
  # try extract seq
  out <- tryCatch(
    seq <- genome[[tmp$seqnames]][tmp$start2:tmp$end2],
    error = function(e){ return(NULL) }
  )
  # check
  if(!is.null(out)){
    # strand
    if(tmp$strand == "-"){
      out <- reverseComplement(out)
    }else{
      out
    }
  }else{
    DNAString()
  }
}) %>% DNAStringSet() -> seq_list

names(seq_list) <- make.unique(genes$gene_id)

seq_list

# filte seq length
fasta_filtered <- seq_list[width(seq_list) > 1]
head(fasta_filtered)
# output fasta format
writeXStringSet(fasta_filtered, filepath = paste0(species, "_3kpromoter.fasta"), format = "fasta")