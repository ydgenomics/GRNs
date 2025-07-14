这个项目是基于时序数据推断调控网络，从分子的调控网络去解释生物的生长发育的事件。
如何去联合分析？

> IReNA环境配置
```shell
source /opt/software/miniconda3/bin/activate
conda create -n IReNA r-base=4.3 -y
conda activate IReNA
conda install bioconda::bioconductor-rsamtools -y
conda install bioconda::bioconductor-chipseeker -y
conda install bioconda::bioconductor-monocle -y
conda install bioconda::bioconductor-rcistarget -y
conda install bioconda::bioconductor-rcy3 -y
conda install bioconda::bioconductor-clusterprofiler -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes -y
conda install conda-forge::r-pbapply -y
conda install conda-forge::r-rocr -y
conda install conda-forge::r-seurat -y # 5.3
conda install bioconda::bioconductor-edger -y
conda install conda-forge::r-furrr -y
conda install conda-forge::r-future -y
Rscript -e 'devtools::install_github("jiang-junyao/IReNA")'
conda install bioconda::bioconductor-annotationdbi -y
conda install bioconda::bioconductor-genie3 -y
conda install conda-forge::r-doparallel -y
conda install conda-forge::r-dorng -y
conda install conda-forge::r-optparse -y
# conda install conda-forge::r-reshape2 -y
```

```R
devtools::install_github("jiang-junyao/IReNA")
library(IReNA)
library(pheatmap)
library(DDRTree)
library(RcisTarget)
library(Rsamtools)
library(GenomicRanges)
library(ChIPseeker)
library(IRanges)
library(Biostrings)
library(igraph)
library(pbapply) # conda install conda-forge::r-pbapply -y
library(RCy3)
library(ROCR) # conda install conda-forge::r-rocr -y
library(Seurat) # conda install conda-forge::r-seurat -y # 5.3
library(ggplot2)
library(dplyr)
library(gridExtra)
library(stats)
library(utils)
library(monocle)
library(VGAM)
library(BiocGenerics)
library(clusterProfiler)
library(edgeR) # conda install bioconda::bioconductor-edger -y
library(limma)
library(rlang)
library(reshape2)
library(stringr)
library(furrr) # conda install conda-forge::r-furrr -y
library(future) # conda install conda-forge::r-future -y
library(purrr)
```

[pp.bulkATAC](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)
```shell
conda create -n htseq -y
conda activate htseq
conda install bioconda::htseq -y

conda create -n samtools -y
conda activate samtools
conda install bioconda::samtools -y

conda create -n bedtools -y
conda activate bedtools
conda install bioconda::bedtools -y

conda create -n macs3 -y
conda activate macs3
conda install -c bioconda macs3 -y
```
[pp.scRNAseq](https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing)

## [bulkATAC+scRNAseq](https://jiang-junyao.github.io/IReNA/scATAC+scRNA)
[samtools]() *Image: chromap*
[MACS3](https://macs3-project.github.io/MACS/index.html) *Image: macs3*
[IReNA]() *Image: IReNA *
[htseq](https://github.com/htseq/htseq) *Image: htseq*
[hint](https://reg-gen.readthedocs.io/en/latest/hint/introduction.html) *Image: rgt-hint*


```shell
conda create -n fastqc -y
conda activate fastqc
conda install bioconda::fastqc -y
conda install bioconda::fastp -y
conda install bioconda::cutadapt -y
conda install bioconda::bowtie2 -y
fastp --version
fastqc --version
cutadapt --version
bowtie2 --version
```

[树棉：*Gossypium arboreum*](https://baike.baidu.com/item/%E6%A0%91%E6%A3%89/1706952?fromModule=search-result_lemma)

计算树棉基因组大小`awk '{sum += length($0)} END {print sum}' /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/genome/genome.fa` **1444625381**

[*Gossypium arboreum* of TF_blinding_motif](https://planttfdb.gao-lab.org/download/motif/Gar_TF_binding_motifs_information.txt)
```R
download.file(
  url = "https://planttfdb.gao-lab.org/download/motif/Gar_TF_binding_motifs_information.txt",
  destfile = "Gar_TF_binding_motifs_information.txt",
  mode = "wb"
)
```

构建cotton的`motif1`
```R
> head(motif1)
  Accession      ID   Name         TFs                             EnsemblID
1    M00001 MYOD_01   MyoD       Myod1                    ENSMUSG00000009471
2    M00002  E47_01    E47        Tcf3                    ENSMUSG00000020167
3    M00004 CMYB_01  c-Myb         Myb                    ENSMUSG00000019982
4    M00005  AP4_01   AP-4       Tfap4                    ENSMUSG00000005718
5    M00006 MEF2_01 MEF-2A Mef2a;Mef2c ENSMUSG00000030557;ENSMUSG00000005583
6    M00007 ELK1_01  Elk-1        Elk1                    ENSMUSG00000009406
```
这个是motif_blinding_tf的信息，最后一列即5列为`gene`
```R
> head(rownames(seurat_object), n=10)
 [1] "Ga14g01907" "Ga05g03180" "Ga13g00005" "Ga14g01554" "Ga14g01555"
 [6] "Ga14g01813" "Ga14g00262" "Ga01g00010" "Ga01g00022" "Ga01g00004"
 ```

 E1_to_Ga的query基因为什么还有不同的后缀，有什么特殊意义吗？？

 rgt-hint需要rgtdata，而一般是人的，需要根据物种进行custom配置[customize-rgt-data-folder](https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html#customize-rgt-data-folder)
 [2019(genome_biology)_Identification of transcription factorbinding sites using ATAC-seq]()

 Motif position weight matirx (PWM)
 [从 gtf 文件构建 orgdb 和 txdb 数据库](https://mp.weixin.qq.com/s/w3FFimm-xF2OY20aoFRcSg)
 *TxDb 用于存储与基因组结构相关的信息，例如染色体、外显子、内含子、启动子、转录本等*
```shell
# txdbmaker只能装在大于等于4.4版本的R里面
conda create -n txdbmaker r-base=4.4 -y
conda activate txdbmaker
conda install conda-forge::r-msigdbr -y
conda install conda-forge::r-tidyverse -y
conda install bioconda::bioconductor-annotationforge -y
conda install bioconda::bioconductor-rtracklayer -y
conda install bioconda::bioconductor-biomart -y
conda install bioconda::bioconductor-txdbmaker -y
 ```

[meme fimo](https://meme-suite.org/meme/doc/install.html)
```shell
# sudo cp /data/work/SCPipelines/meme-5.5.8.tar.gz .
wget https://meme-suite.org/meme/meme-software/5.5.8/meme-5.5.8.tar.gz
cd ~
tar zxf /data/work/SCPipelines/meme-5.5.8.tar.gz
cd meme-5.5.8
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install

export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.8:$PATH
export PATH=$PATH:~/meme/bin
```

[【FDU生统实践课-3】单细胞ATAC-seq全分析流程 30分钟入门](https://www.bilibili.com/video/BV1Be4y127Ho/?spm_id_from=333.337.search-card.all.click&vd_source=5600c17ea3ce6334fe6d9c0d3cd99627)