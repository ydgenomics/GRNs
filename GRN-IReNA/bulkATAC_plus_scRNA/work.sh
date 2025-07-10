mkdir test
cd test
mkdir fimo
mkdir ./fimo/output
work_path=$(pwd)
echo "Currently working path: $work_path"

time_rds="/data/input/Files/ResultData/Workflow/W202506200040331/result/seu_day1.rds_seurat_with_time.rds"
input_bams_txt="/data/work/test/input_bams.txt"
# bam_inputs <- scan(input_bams_txt, what = character(), sep = ",")
tf2motif_txt="/data/work/SCPipelines/pp_scRNA/tf_blinding_motif_ga.txt"
genie3_rdata="/data/input/Files/ResultData/Workflow/W202507020005242/regulatory_relationships.RData"
# footprints_bed="/data/work/SCPipelines/pp_bulkATAC/10.footprints/filtered_footprints0.bed" 
footprints_bed="/data/work/SCPipelines/pp_bulkATAC/10.footprints/filtered_footprints_0.6.bed"
genome_fa="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa"
pwm_dir="/data/work/SCPipelines/Gar_individual_motif2/"
peaks_bed="/data/work/test0707/08_differential_peaks/peaks.bed"
txdb_sqlite="/data/work/SCPipelines/gtf2txdb/txdb.sqlite"
get_merged_fasta_R='/data/work/SCPipelines/all/get_merged_fasta.R'
annotation_R="/data/work/SCPipelines/all/annotation.R"

/opt/software/miniconda3/envs/IReNA/bin/Rscript /data/work/test/bulkATAC_scRNA.R \
$time_rds $input_bams_txt $tf2motif_txt $genie3_rdata $footprints_bed $genome_fa \
$pwm_dir $peaks_bed $txdb_sqlite $work_path $get_merged_fasta_R $annotation_R
