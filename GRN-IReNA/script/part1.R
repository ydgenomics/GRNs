# Part 1: Analyze scRNA-seq or bulk RNA-seq data to get basic regulatory relationships
###Load the test data
seurat_with_time <- readRDS('seurat_with_time.rds')
###Get expression profiles ordered by pseudotime
expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 50)
###Filter noise and logFC in expression profile
expression_profile_filter <- filter_expression_profile(expression_profile, FC=0.01)
###K-means clustering
clustering <- clustering_Kmeans(expression_profile_filter, K1=4)

clustering[1:5,1:5]

pdf('plot_kmeans_pheatmap.pdf', width=10, height=10)
plot_kmeans_pheatmap(clustering,ModuleColor1 = c('#67C7C1','#5BA6DA','#FFBF0F','#C067A9'))
dev.off()

# GENIE3[optional] to infer regulatory relationships
weightMat <- GENIE3(as.matrix(seurat_with_time@assays$RNA@data),nCores = 50)
weightMat <- getLinkList(weightMat)
regulation <- weightMat[weightMat[,3]>0.0002,]
### add regulation type for each gene pair
regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS,regulation)
### check whether source genes are transcription factors
motifTF <- c()
for (i in 1:nrow(motif1)) {
  TF <- strsplit(motif1[i,5],';')[[1]]
  motifTF <- c(motifTF,TF)
}
regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]

# Person’s correlation[optional] to filter regulatory relationships
regulatory_relationships <- get_cor(Kmeans_clustering_ENS, motif = motif1, correlation_filter = 0.6, start_column = 4)

# Part 2: Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)

###merge footprints whose distance is less than 4
filtered_footprints <- read.table('footprints.bed',sep = '\t')
fastadir <- 'Genome/hg38.fa' 
merged_fasta <- get_merged_fasta(filtered_footprints,fastadir)
write.table(merged_fasta,'merged_footprints.fasta',row.names=F,quote=F,col.names=F)

### Identify differentially expressed genes related motifs
motif1 <- motifs_select(Tranfac201803_Hs_MotifTFsF, rownames(Kmeans_clustering_ENS)) ###Kmeans_clustering_ENS was obtained in part1
### run find_motifs()
fimodir <- 'fimo'
outputdir1 <- '/public/home/user/fimo/output/'
outputdir <- '/public/home/user/fimo/output/'
motifdir <- '/public/home/user/fimo/Mememotif/'
sequencedir <- '/public/home/user/fimo/merged_footprints.fasta'
find_motifs(motif1,step=20,fimodir, outputdir1, outputdir, motifdir, sequencedir)
### run fimo_all script in shell
shell_code <- paste0('sh ',outputdir1,'Fimo_All.sh')
system(shell_code,wait=TRUE)
### delete shell scripts
shell_code2 <- paste0('rm ',outputdir1,'Fimo*.sh')
system(shell_code2,wait=TRUE)

###Combine all footprints of motifs
combined <- combine_footprints(outputdir)
peaks <- read.delim('differential_peaks.bed')
overlapped <- overlap_footprints_peaks(combined,peaks)

###get footprint-related genes
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # Build species specify TxDb
list1 <- get_related_genes(overlapped,txdb = txdb,motif=Tranfac201803_Hs_MotifTFsF,Species = 'Hs')
###Get candidate genes/TFs-related peaks
list2 <- get_related_peaks(list1,Kmeans_clustering_ENS)
### output filtered footprints
write.table(list2[[1]],'filtered_footprints.bed', quote = F, row.names = F, col.names = F, sep = '\t')

### run samtools in shell[optional]
shell_code1 <- 'samtools view -hb -L filtered_footprint.bed SSC_patient1.bam > SSC1_filter.bam'
shell_code2 <- 'samtools view -hb -L filtered_footprint.bed SSC_patient2.bam > SSC2_filter.bam'
shell_code3 <- 'samtools view -hb -L filtered_footprint.bed esc.bam > esc_filter.bam'
system(shell_code1,,wait=TRUE)
system(shell_code2,,wait=TRUE)
system(shell_code3,,wait=TRUE)

# [optional]
library(ATACseqQC)
library(Rsamtools)
bamfilepath1 <- 'SSC1_filter.bam'
bamfilepath2 <- 'SSC2_filter.bam'
bamfilepath3 <- 'esc_filter.bam'
indexBam(bamfilepath1)
gal1 <- readBamFile(bamfilepath1, tag=tags, asMates=TRUE, bigFile=TRUE)
gal2 <- readBamFile(bamfilepath2, tag=tags, asMates=TRUE, bigFile=TRUE)
gal3 <- readBamFile(bamfilepath3, tag=tags, asMates=TRUE, bigFile=TRUE)
galout1 <- shiftGAlignmentsList(gal, 'SSC1_filter_shift.bam')
galout2 <- shiftGAlignmentsList(gal, 'SSC2_filter_shift.bam')
galout3 <- shiftGAlignmentsList(gal, 'esc_filter_shift.bam')

### calculate cuts of each each position in footprints
bamfilepath1 <- 'SSC1_filter.bam'
bamfilepath2 <- 'SSC2_filter.bam'
bamfilepath3 <- 'esc_filter.bam'
### set parameter 'workers' to make this function run in parallel
cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1,bedfile = list2[[1]],workers = 40,index_bam = T)
cuts2 <- cal_footprint_cuts(bamfilepath = bamfilepath2,bedfile = list2[[1]],workers = 40,index_bam = T)
cuts3 <- cal_footprint_cuts(bamfilepath = bamfilepath3,bedfile = list2[[1]],workers = 40,index_bam = T)
cut_list <- list(cuts1,cuts2,cuts3)
### get related genes of footprints with high FOS
potential_regulation <- Footprints_FOS(cut_list,list2[[2]], FOS_threshold = 0.1)
### Use information of footprints with high FOS to refine regulatory relationships
filtered_regulatory <- filter_ATAC(potential_regulation,regulatory_relationships)

# Part 3: Regulatory network analysis and visualization
# filtered_regulatory_relationships Kmeans_clustering_ENS
TFs_list <- network_analysis(filtered_regulatory_relationships,Kmeans_clustering_ENS,TFFDR1 = 10,TFFDR2 = 10)
