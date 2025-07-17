# Based on [IReNA](https://github.com/jiang-junyao/IReNA) do GRN analysis
**Brief**:  There are several cell states involved in cell development or disease occurrence (e.g., progenitor, precursor, immature, and mature), each state maintained by a unique gene program (gene regulatory modules). Decoding the inter- or intra-regulatory mechanisms among these modules can further elucidate the key mechanisms that regulate cell state transitions, including identifying key transcription factors that regulate cell fate decisions or cell differentiation. Most current gene regulatory network (GRN) analysis methods focus on intra-modular regulations; they select all cell states or single cell states to construct GRNs and neglect inter-modular regulations.

IReNA can address this gap by identifying transcription factors (TFs) that regulate other modules and inferring inter-modular interactions through hypergeometric tests. For instance, if IReNA identifies TF A from module a significantly activating module b, we can infer that TF A may regulate the differentiation of the progenitor state into the precursor state. In a second case, if IReNA identifies TF B from module c significantly repressing module d, we can infer that TF B represses the differentiation process from the immature state to the mature state.

# Overview
![workflow of IReNA](https://github.com/jiang-junyao/IReNA/blob/master/docs/Readme%20figure/workflow_new.jpg)

- scRNA: Regulatory network analysis through only scRNA-seq data [IReNA](https://jiang-junyao.github.io/IReNA/only-scRNA) [ydgenomics]()
- scRNA+scATAC: [Regulatory network analysis through intergrating scRNA-seq data and scATAC-seq data](https://jiang-junyao.github.io/IReNA/scATAC+scRNA) [ydgenomics]()
- scRNA+bulkATAC: [Regulatory network analysis through intergrating scRNA-seq data and bulk ATAC-seq data](https://jiang-junyao.github.io/IReNA/bulk-ATAC+scRNA) [ydgenomics]()

# Input of Preparation
- Requiried
  1. Genome(.fa .fai)
  2. Annotation of genome(.gtf)
  3. TF_blinding_motif(downloaded from [PlantTFDB]())
  4. PWM(position weight matrix)(downloaded from PlantTFDB)
  5. scRNA(.rds)

- Optional
  1. scRNA+bulkATAC Requiried: rgtdata()
  2. scRNA+scRNA

# Scripts & Piplines
ATAC-seq preprocessing pipeline
  - Step 1: convert sra file to fastq file `sratools` `.fastq`
  - Step 2: Quality control `fastqc, fastp, cutadapt` `.fastq`
    - `fastqc` to check reads quality
    - `fastp` is a user friendly software which can remove adapator and low quality reads
    - `cutadapt` to remove adaptor
    - > what is 'adapator sequence'?
  - Step 3: Mapping `bowtie2` `.sam`
    - `bowtie2-build` create bowtie2 index with reference genome
    - `bowtie2` get result `.sam` of mapping
  - **Step 4**: Filter low-quality reads (MAPQ score < 10), then remove duplicate with `samtools` `.bam`
  - Step 5: Call peaks with `macs2` `.txt`
  - Step 6: Merge all peaks and generate gtf file `library(IReNA)` `.gtf`
  - Step 7: Calculate raw counts withn htseq `htseq-count` `.txt`
  - **Step 8**: Combine all samples, and calculate differential peaks according to edgeR `library(IReNA)` `differential_peaks.bed, peaks.bed`
  - Step 9: Merge bam files `samtools`
  - **Step 10**: Calculate footprints of merged bam files with HINT `rgt-hint` `differential_peaks.bed, filtered_footprints.bed`

# Output & Plot
