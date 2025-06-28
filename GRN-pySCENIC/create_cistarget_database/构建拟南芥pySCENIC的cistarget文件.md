# 构建拟南芥pySCENIC的cistarget文件

**Reference:** [中文](https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ) [github](https://github.com/aertslab/create_cisTarget_databases)
**最终的目录内容**：
```shell
tree -L 2
.
├── Athaliana_TAIR10.54.dna.fa
├── Athaliana_TAIR10.54.gtf
├── Ath_pep.fas
├── Ath_TF_binding_motifs_information.txt
├── Ath_TF_binding_motifs.meme
├── Ath_TF_list.txt
├── output
│   ├── AT_3kpromoter.fasta
│   ├── Ath_TF_list.txt
│   ├── AT_motif_dir
│   ├── AT_motifs_id.txt
│   ├── AT.motifs_vs_regions.scores.feather
│   ├── AT.regions_vs_motifs.rankings.feather
│   ├── AT.regions_vs_motifs.scores.feather
│   └── AT_TF_binding_motifs_information.tbl
└── script
    ├── create_cistarget_motif_databases_yd.py
    ├── extract_promoters.R
    └── get_AT.file.sh
```
output里面用于下游分析的三个文件：`Ath_TF_list.txt`, `AT_TF_binding_motifs_information.tbl`, `AT.regions_vs_motifs.rankings.feather`。
[百度网盘下载，提取码：1111](https://pan.baidu.com/s/15KNB1LjBhSqk8VvrD27AdA)

## 1.下载拟南芥的TF列表和对应的蛋白质序列文件,TF_motif对应信息和motif.meme文件
![AT_gene_protein](../png/AT_gene_name.png)

[**Download:** Arabidopsis thaliana TF list](https://planttfdb.gao-lab.org/download/TF_list/Ahy_TF_list.txt.gz)
[**Download:** Arabidopsis thaliana Protein sequences](https://planttfdb.gao-lab.org/download/seq/Ath_pep.fas.gz)

![AT_tfmotif_meme](../png/AT_tfmotif_meme.png)

[**Download:** Arabidopsis thaliana TF_binding_motifs](https://planttfdb.gao-lab.org/download/motif/Ath_TF_binding_motifs_information.txt)
[**Download:** Arabidopsis thaliana motifs.meme](https://planttfdb.gao-lab.org/download/motif/Ath_TF_binding_motifs.meme.gz)
```R
# Suggestion download Arabidopsis thaliana TF_binding_motifs uses your PC
download.file(
  url = "https://planttfdb.gao-lab.org/download/motif/Ath_TF_binding_motifs_information.txt",
  destfile = "Ath_TF_binding_motifs_information.txt",
  mode = "wb"
)
```

## 2.处理拟南芥的TF_motif和motif.meme文件 
[get_AT.file.sh](../create_cistarget_database/get_AT.file.sh)
>AT_motif_dir
>AT_motifs_id.txt
>AT_TF_binding_motifs_information.tbl


## 3.获取拟南芥基因组的启动子序列，并根据motif信息构建cistarget的database
[extract_promoters.R](../create_cistarget_database/extract_promoters.R)拿到`AT_3kpromoter.fasta`
```sehll
Rscript extract_promoters.R \
--gtf /data/work/0.peanut/GRN/AT/Athaliana_TAIR10.54.gtf \
--fasta /data/work/0.peanut/GRN/AT/Athaliana_TAIR10.54.dna.fa \
--species AT
```

使用[create_cistarget_motif_databases_yd.py](../create_cistarget_database/create_cistarget_motif_databases_yd.py)构建cistarget_database，其中`AT.regions_vs_motifs.rankings.feather`用做后面分析
```shell
python /script/grn_pyscenic/create_cistarget_motif_databases_yd.py \
    -f /data/work/0.peanut/GRN/AT/output/AT_3kpromoter.fasta \
    -M /data/work/0.peanut/GRN/AT/output/AT_motif_dir/ \
    -m /data/work/0.peanut/GRN/AT/output/AT_motifs_id.txt \
    -t 8 \
    -o AT
```

## 保存没有后缀以AT开头的名字作为AT_TF文件，每一个基因占一行
```R
file_path <- "/data/work/0.peanut/GRN/AT/Ath_TF_list.txt"
tf_list <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
head(tf_list)
gene_ids <- tf_list$Gene_ID
# 写入新的文本文件，每个值占一行
output_file <- "Ath_TF_list.txt"
write(gene_ids, file = output_file, sep = "\n")
```