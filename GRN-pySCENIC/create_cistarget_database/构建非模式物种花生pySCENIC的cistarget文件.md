# 构建非模式物种花生pySCENIC的cistarget文件

**Reference:** [中文](https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ) [github](https://github.com/aertslab/create_cisTarget_databases)

**最终的目录内容**：
```shell

```

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
>AT_TF_binding_motifs_information.tbl.txt


## 3.准备花生的基因组相关文件用于比对和建cistarget库
```shell
genome_fasta="/data/input/Files/husasa/Ref/arahy.Tifrunner.gnm2.J5K5.genome_main.fa"
genome_gtf="/data/work/0.peanut/GRN/output/updated_gtf_file_standard.gtf"
protein_fasta="/data/work/0.peanut/GRN/input/arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa"
head -n 3 $genome_fasta
#>arahy.Tifrunner.gnm2.chr01
#AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAACCTAAACCCTAAACCCTAAACCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAC
#CCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCC
head -n 3 $genome_gtf
#arahy.Tifrunner.gnm2.chr01      Mikado_loci     transcript      19126   25719   16      +       NA      transcript_id "arahy.Tifrunner.gnm2.ann2.Ah01g000200.1"; gene_id "arahy.Tifrunner.gnm2.ann2.Ah01g000200.1";
#arahy.Tifrunner.gnm2.chr01      Mikado_loci     exon    19126   19514   NA      +       NA      transcript_id "arahy.Tifrunner.gnm2.ann2.Ah01g000200.1"; gene_id "arahy.Tifrunner.gnm2.ann2.Ah01g000200.1";
#arahy.Tifrunner.gnm2.chr01      Mikado_loci     exon    20456   20598   NA      +       NA      transcript_id "arahy.Tifrunner.gnm2.ann2.Ah01g000200.1"; gene_id "arahy.Tifrunner.gnm2.ann2.Ah01g000200.1";
head -n 3 $protein_fasta
#>arahy.Tifrunner.gnm2.ann2.Ah01g000200.1
#MSVAADSPIHSSSSDDFIAYLDDALAASSPDASSDKEVENQDELESGRIKRCKFESAEETEESTSEGIVKQNLEEYVCTHPGSFGDMCIRCGQKLDGESGVTFGYIHKGLRLHDEEISRLRNTDVKNLLIRKKLYLILDLDHTLLNSTHLAHLNSEELHLISQADSLGDVSKGSLFKLDKMHMMTKLRPFVRTFLKEASEMFEMYIYTMGDRPYALEMAKLLDPLGEYFNAKVISRDDGTQKHQKGLDIVLGQESAVVILDDTEHAWVKHKDNLILMERYHFFGSSCRQFGFNCKSLAELKSDEDEAEGALTKILKVLKQVHSKFFDELKEDIAERDVRQVLKSVRREVLSGCVVVFSRIFHGALPPLRQMAEQLGATCLMELDPSVTHVVATDAGTEKARWAVKEKKFLVHPRWIEAANYFWEKQPEENFVLKKKQ
#>arahy.Tifrunner.gnm2.ann2.Ah01g000400.1
```

## 4.获取花生基因组的启动子序列，并根据motif信息构建cistarget的database
[extra_promoters.R](extract_promoters.R)拿到`3kpromoter.fasta`
使用[create_cistarget_motif_databases_yd.py](./create_cistarget_database/create_cistarget_motif_databases.py)构建cistarget_database，其中`peanut.regions_vs_motifs.rankings.feather`用做后面分析
```shell
python /data/work/0.peanut/GRN/create_cistarget_motif_databases_yd.py \
-f /data/work/0.peanut/GRN/output/3kpromoter.fasta \
-M /data/work/0.peanut/GRN/output/motif_dir/ \
-m /data/work/0.peanut/GRN/output/motifs_id.txt \
-t 30 \
-o peanut
```

## 5.安装blastp并做两个物种的蛋白质比对找[同源基因](../../Ortho_gene/)
**Reference：** [史上最详细的blast安装附视频](https://mp.weixin.qq.com/s/rEBqjN-fGOp_loTmyEuMJA)
```shell
########## Install blastp ##########
cd /software
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
# environment location: /software/ncbi-blast-2.16.0+/bin
vim ~/.bashrc
export PATH=/software/ncbi-blast-2.16.0+/bin:$PATH
source ~/.bashrc
blastp -h
rm ncbi-blast-2.16.0+-x64-linux.tar.gz
```
```shell
########## Using blastp ##########
subject_fasta="/data/work/0.peanut/GRN/input/Ath_pep.fas"
query_fasta="/data/work/0.peanut/GRN/input/arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa"
#cd /data/work/0.peanut/GRN/input/arabidopsis_db
# Make database of AT/Others
makeblastdb -in $subject_fasta -dbtype prot -out arabidopsis_db
# Query input fasta
blastp -query $query_fasta -db arabidopsis_db -out blastp_results.txt -outfmt 6 -evalue 1e-5
```

## 6.处理比对信息拿到花生和拟南芥一对一的关系文件，修改拟南芥的tbl和TF_list，拿到peanut的tbl和TF_list
[get_one2one_map2AT.file.py](./create_cistarget_database/get_one2one_map2AT.file.py)

## 如果基因组基因名和单细胞数据矩阵的基因名不一致，还需要读取gtf文件和矩阵基因名处理，可参考[T-peanut/change_gname_gtf](https://github.com/ydgenomics/T-peanut/blob/main/change_gname_gtf.R)