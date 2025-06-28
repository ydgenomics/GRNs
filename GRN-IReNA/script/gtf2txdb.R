# Date: 20250617 # Title: gtf2txdb.R
# Ref:  [从 gtf 文件构建 orgdb 和 txdb 数据库](https://mp.weixin.qq.com/s/w3FFimm-xF2OY20aoFRcSg)
# /software/miniconda/envs/txdbmaker/bin/R

library(msigdbr)
library(tidyverse)
library(AnnotationForge)
library(AnnotationDbi)
library(rtracklayer)
library(biomaRt)
library(txdbmaker)

# make txdb from gtf https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromGFF.html
# [Taxonomy: Gossypium arboreum](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=29729)
txdb <- makeTxDbFromGFF(file = "/data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/ga.gtf", 
                        organism = "Gossypium arboreum")
# Import genomic features from the file as a GRanges object ... OK
# Prepare the 'metadata' data frame ... OK
# Make the TxDb object ... OK

txdb
# TxDb object:
# # Db type: TxDb
# # Supporting package: GenomicFeatures
# # Data source: /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/SCPipelines/pp_bulkATAC/build_rgtdata/ga/ga.gtf
# # Organism: Gossypium arboreum
# # Taxonomy ID: 29729
# # miRBase build ID: NA
# # Genome: NA
# # Nb of transcripts: 94545
# # Db created by: txdbmaker package from Bioconductor
# # Creation time: 2025-06-17 12:22:28 +0800 (Tue, 17 Jun 2025)
# # txdbmaker version at creation time: 1.2.0
# # RSQLite version at creation time: 2.4.1
# # DBSCHEMAVERSION: 1.2

# Save txdb as .sqlite
saveDb(txdb, file = "/data/work/test/txdb.sqlite")
# txdb <- loadDb("/data/work/test/txdb.sqlite")