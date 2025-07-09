################################## 10_footprints ##################################
mkdir 10_footprints
cd 10_footprints
bam="../09_merge_bams/merged_all.bam"
peaks_bed="../08_differential_peaks/peaks.bed"

# # Building rgtdata for TAIR10(Maybe another organism)
# genome="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa"
# genome_fai="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa.fai"
# chromosome_sizes="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/chrom.sizes.ga"
# gene_regions="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genes_ga.bed"
# annotation="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/ga.gtf"
# gene_alias="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/alias_ga.txt"

# sudo mkdir ~/rgtdata/tair10
# sudo cp $genome ~/rgtdata/tair10/genome_tair10_ensembl_release_51.fa
# sudo cp $genome_fai ~/rgtdata/tair10/genome_tair10_ensembl_release_51.fa.fai
# sudo cp $chromosome_sizes ~/rgtdata/tair10/chrom.sizes.tair10
# sudo cp $gene_regions ~/rgtdata/tair10/genes_tair10.bed
# sudo cp $annotation ~/rgtdata/tair10/Arabidopsis_thaliana.TAIR10.51.gtf
# sudo cp $gene_alias ~/rgtdata/tair10/alias_tair10.txt

# rgt-hint footprinting --atac-seq --paired-end --organism=tair10 \
# ../09_merge_bams/merged_all.bam ../08_differential_peaks/differential_peaks.bed

rgt-hint footprinting --atac-seq --paired-end --organism=tair10 \
$bam $peaks_bed

R_CODE='
### R code
# footprints <- read.table("footprints.bed",sep="\t",header = T)
footprints <- read.table("footprints.bed", sep="\t", header = FALSE)
footprints_80th <- footprints[footprints$V5 > quantile(footprints$V5, 0.8),] # `quantile`
footprints_80th <- footprints_80th[,c(1,2,3,5)]
write.table(footprints_80th, "filtered_footprints.bed", quote=F, sep = "\t", row.names = F, col.names = F)
footprints_80th <- footprints[footprints$V5 > quantile(footprints$V5, 0.6),] # `quantile`
footprints_80th <- footprints_80th[,c(1,2,3,5)]
write.table(footprints_80th, "filtered_footprints_0.6.bed", quote=F, sep = "\t", row.names = F, col.names = F)
'

echo "$R_CODE" | /opt/software/miniconda3/envs/IReNA/bin/R --vanilla