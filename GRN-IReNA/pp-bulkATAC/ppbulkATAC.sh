# bams_txt="/data/work/test0629/bams.txt"
# names_txt="/data/work/test0629/names.txt"
# genome_size=1444625381

# sh /data/work/test0710/ppbulkATAC.sh bams.txt names.txt 1444625381
bams_txt=$1
names_txt=$2
genome_size=$3

# bams_txt="/data/work/test0710/bams.txt"
# names_txt="/data/work/test0710/names.txt"
# genome_size=1444625381

# 读取文件内容并按逗号分割
IFS=',' read -r -a names <<< "$(cat "$names_txt")"
IFS=',' read -r -a bams <<< "$(cat "$bams_txt")"

# 获取数组长度
len_names=${#names[@]}
len_bams=${#bams[@]}

# 确保两个数组长度一致
if [ "$len_names" -ne "$len_bams" ]; then
  echo "Error: The number of names and outputs does not match."
  exit 1
fi

# samtools view /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-2_bwa_rmdup.bam | head -n 10
################################## 05_call_peaks ##################################
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "05_call_peaks" >> log.txt
echo "Start-time of 05_call_peaks: $start_time" >> log.txt

mkdir 05_call_peaks
cd 05_call_peaks

# 遍历数组并执行 macs3 callpeak 命令
for ((i=0; i<${#names[@]}; i++)); do
    name=${names[i]}
    bam=${bams[i]}
    log_file="macs3_log_${name}.txt"
    /opt/software/miniconda3/envs/macs3/bin/macs3 callpeak -t "$bam" --nomodel --shift -100 --extsize 200 -g $genome_size -n "$name" > "$log_file"
done

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End-time of 06_merge_peaks: $end_time" >> ../log.txt
# Duration time
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$((end_seconds - start_seconds))
echo "Duration seconds of 05_call_peaks: $duration" >> ../log.txt
################################## 06_merge_peaks ##################################
# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "06_merge_peaks" >> ../log.txt
echo "Start-time of 06_merge_peaks: $start_time" >> ../log.txt

mkdir ../06_merge_peaks
cd ../06_merge_peaks

# 定义 R 代码
R_CODE='
library(IReNA)
# 读取文件内容
peak_files <- "../names.txt"
peak_files <- readLines(peak_files)
peak_files <- unlist(strsplit(peak_files, ","))

peak_files <- paste0("../05_call_peaks/", peak_files, "_peaks.narrowPeak")
peak_list <- lapply(peak_files, function(f) read.delim(f, header = FALSE))
all_peak <- do.call(rbind, peak_list)
peaks_merged_gtf <- generate_peak_gtf(all_peak)

head(peaks_merged_gtf)

write.table(peaks_merged_gtf, "peaks_merged.gtf", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
'

# 运行 R 代码
echo "$R_CODE" | /opt/software/miniconda3/envs/IReNA/bin/R --vanilla

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End-time of 06_merge_peaks: $end_time" >> ../log.txt
# Duration time
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$((end_seconds - start_seconds))
echo "Duration seconds of 06_merge_peaks: $duration" >> ../log.txt
################################## 07_calculate_counts ##################################
# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "07_calculate_counts" >> ../log.txt
echo "Start-time of 07_calculate_counts: $start_time" >> ../log.txt

mkdir ../07_calculate_counts
cd ../07_calculate_counts
# 遍历数组并执行 macs3 callpeak 命令
for ((i=0; i<${#names[@]}; i++)); do
    name=${names[i]}
    bam=${bams[i]}
    counts_file="${name}_Counts.txt"
    htseq-count -r pos --idattr gene_id --stranded no -a 10 -f bam "$bam" ../06_merge_peaks/peaks_merged.gtf > $counts_file
done

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End-time of 07_calculate_counts: $end_time" >> ../log.txt
# Duration time
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$((end_seconds - start_seconds))
echo "Duration seconds of 07_calculate_counts: $duration" >> ../log.txt
################################## 08_differential_peaks ##################################
# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "08_differential_peaks" >> ../log.txt
echo "Start-time of 08_differential_peaks: $start_time" >> ../log.txt

mkdir ../08_differential_peaks
cd ../08_differential_peaks

R_CODE='
library(IReNA)
count_files <- "../names.txt"
count_files <- readLines(count_files)
count_files <- unlist(strsplit(count_files, ","))
count_files <- paste0("../07_calculate_counts/", count_files, "_Counts.txt")
count_list <- list()
for (file in count_files) {
    count_list[[file]] <- read.delim(file, header = FALSE)
}
head(count_list[[1]])
count_all <- count_list[[1]]
if (length(count_list) > 1){
    for (i in 2:length(count_list)) {
        count_all <- cbind(count_all, count_list[[i]][, 2])
    }
}
head(count_all)
peaks_gtf <- read.delim("../06_merge_peaks/peaks_merged.gtf", header=FALSE)
# added by yd
merge_sort_count <- function (count_all, gtf) {
    # validInput(count, "motif", "df")
    # validInput(gtf, "gtf", "df")
    count_all <- count_all[1:(nrow(count_all) - 5), ]
    gtf$num <- 1:nrow(gtf)
    gtf <- gtf[order(gtf[, 10]), ]
    count_all <- cbind(count_all, gtf[, c(1, 4, 5, 11)])
    count_all <- count_all[order(count_all[, ncol(count_all)]), 
        ]
    rownames(count_all) <- count_all[, 1]
    count_all <- count_all[, -1]
    count_all <- count_all[, c((ncol(count_all) - 3):(ncol(count_all) - 
        1), 1:(ncol(count_all) - 4))]
    return(count_all)
}
merged_count <- merge_sort_count(count_all, peaks_gtf)
head(merged_count) ### first column is chromosome, second column is peak start site, third column is peak end site, the fourth through sixth columns are counts
# merged_count saved as .csv file
write.csv(merged_count, file = "merged_count.csv", row.names = FALSE)
write.table(merged_count[, c(1,2,3)], "peaks.bed", quote = F, row.name = F, col.names = F, sep = "\t")
'
# 运行 R 代码
echo "$R_CODE" | /opt/software/miniconda3/envs/IReNA/bin/R --vanilla

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End-time of 08_differential_peaks: $end_time" >> ../log.txt
# Duration time
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$((end_seconds - start_seconds))
echo "Duration seconds of 08_differential_peaks: $duration" >> ../log.txt
################################## 09_merge_bams ##################################
# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "09_merge_bams" >> ../log.txt
echo "Start-time of 09_merge_bams: $start_time" >> ../log.txt

mkdir ../09_merge_bams
cd ../09_merge_bams

/opt/software/miniconda3/envs/samtools/bin/samtools merge merged_all.bam "${bams[@]}"

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End-time of 09_merge_bams: $end_time" >> ../log.txt
# Duration time
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$((end_seconds - start_seconds))
echo "Duration seconds of 09_merge_bams: $duration" >> ../log.txt
################################## 10_footprints ##################################
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "10_footprints" >> ../log.txt
echo "Start-time of 10_footprints: $start_time" >> ../log.txt

mkdir ../10_footprints
cd ../10_footprints

# # genome="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa"
# # genome_fai="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genome_ga.fa.fai"
# # chromosome_sizes="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/chrom.sizes.ga"
# # gene_regions="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/genes_ga.bed"
# # annotation="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/ga.gtf"
# # gene_alias="/data/work/SCPipelines/pp_bulkATAC/build_rgtdata/ga/alias_ga.txt"

# # sudo mkdir ~/rgtdata/tair10
# # sudo cp $genome ~/rgtdata/tair10/genome_tair10_ensembl_release_51.fa
# # sudo cp $genome_fai ~/rgtdata/tair10/genome_tair10_ensembl_release_51.fa.fai
# # sudo cp $chromosome_sizes ~/rgtdata/tair10/chrom.sizes.tair10
# # sudo cp $gene_regions ~/rgtdata/tair10/genes_tair10.bed
# # sudo cp $annotation ~/rgtdata/tair10/Arabidopsis_thaliana.TAIR10.51.gtf
# # sudo cp $gene_alias ~/rgtdata/tair10/alias_tair10.txt

# # rgt-hint footprinting --atac-seq --paired-end --organism=tair10 \
# # ../09_merge_bams/merged_all.bam ../08_differential_peaks/differential_peaks.bed

rgt-hint footprinting --atac-seq --paired-end --organism=tair10 \
../09_merge_bams/merged_all.bam ../08_differential_peaks/peaks.bed

# # rgt-hint footprinting --atac-seq --paired-end --organism=tair10 \
# # /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/test0707/09_merge_bams/merged_all.bam \
# # /data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/test0707/08_differential_peaks/peaks.bed

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
# library(IReNA)
# setwd("/data/work/test0707/10_footprints")
# footprints <- read.table("footprints.bed", sep="\t", header = FALSE)
# footprints_80th <- footprints[footprints$V5 > quantile(footprints$V5, 0.9),] # `quantile`
# footprints_80th <- footprints_80th[,c(1,2,3,5)]
# write.table(footprints_80th, "filtered_footprints_0.9.bed", quote=F, sep = "\t", row.names = F, col.names = F)

echo "$R_CODE" | /opt/software/miniconda3/envs/IReNA/bin/R --vanilla

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End-time of 10_footprints: $end_time" >> ../log.txt
# Duration time
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$((end_seconds - start_seconds))
echo "Duration seconds of 10_footprints: $duration" >> ../log.txt