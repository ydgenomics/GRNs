
# bams_txt="/data/work/test0629/bams.txt"
# names_txt="/data/work/test0629/names.txt"
# genome_size=1444625381

bams_txt=$1
names_txt=$2
genome_size=$3

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
mkdir 05_call_peaks
cd 05_call_peaks

# 遍历数组并执行 macs3 callpeak 命令
for ((i=0; i<${#names[@]}; i++)); do
    name=${names[i]}
    bam=${bams[i]}
    log_file="macs3_log_${name}.txt"
    /opt/software/miniconda3/envs/macs3/bin/macs3 callpeak -t "$bam" --nomodel --shift -100 --extsize 200 -g $genome_size -n "$name" > "$log_file"
done

################################## 06_merge_peaks ##################################
mkdir ../06_merge_peaks
cd ../06_merge_peaks

# 定义 R 代码
R_CODE='
library(IReNA)
# 读取文件内容
peak_files <- "/data/work/test0629/names.txt"
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

################################## 07_calculate_counts ##################################
mkdir ../07_calculate_counts
cd ../07_calculate_counts
# 遍历数组并执行 macs3 callpeak 命令
for ((i=0; i<${#names[@]}; i++)); do
    name=${names[i]}
    bam=${bams[i]}
    counts_file="${name}_Counts.txt"
    htseq-count -r pos --idattr gene_id --stranded no -a 10 -f bam "$bam" ../06_merge_peaks/peaks_merged.gtf > $counts_file
done

################################## 08_differential_peaks ##################################
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
for (i in 2:length(count_list)) {
    count_all <- cbind(count_all, count_list[[i]][, 2])
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
write.table(merged_count[, c(1,2,3)], "peaks.bed", quote = F, row.name = F, col.names = F, sep = "\t")
# # group1 <- c(1,1,2)  ### set group for three sample
# group1 <- c(1,1,1,2) # added by yd
# # differential_peaks1 <- diff_peaks(merged_count[,4:6], group1) ### identify differential peaks
# differential_peaks1 <- diff_peaks(merged_count[,4:7], group1) # added by yd
# head(differential_peaks1)
# differential_peaks2 <- merged_count[differential_peaks1[,3]<0.05,] ### filter peaks whose FDR is more than 0.05
# write.table(differential_peaks2[, c(1,2,3)], "differential_peaks.bed", quote = F, row.name = F, col.names = F, sep = "\t")

### If you only have one sample, run the following codes
# SSC_patient1_Counts <- read.delim("SSC_patient1_Counts.txt", header=FALSE)
# peaks_gtf <- read.delim("peak_merged.gtf", header=FALSE)
# merged_count <- merge_sort_count(SSC_patient1_Counts, peaks_gtf)
# write.table(merged_count[, c(1,2,3)], "peaks.bed", quote = F, row.name = F, col.names = F, sep = "\t")
'
# 运行 R 代码
echo "$R_CODE" | /opt/software/miniconda3/envs/IReNA/bin/R --vanilla
################################## 09_merge_bams ##################################
mkdir ../09_merge_bams
cd ../09_merge_bams

/opt/software/miniconda3/envs/samtools/bin/samtools merge merged_all.bam "${bams[@]}"