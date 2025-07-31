# Using [Gene2Role](https://github.com/liushu2019/Gene2Role) reveals the different topolgy genes(DTGs) [GRN-gene2role](https://github.com/ydgenomics/GRNs/tree/main/GRN-gene2role)
- **Log**
  - 250731 完成第一版流程搭建，后面补充DTGs的模块分析和富集分析


# Input
- rds RNA有counts
- txt 文本文件，三列分别对应 `from` `to` `orientation` orientation列只包含1和-1，无行名，无列名
- cluster_key 分群信息对应的键
- cluster_value 取需要对比的群，用`,`连接；如果全部都关注则输入`all`
- n_hvg 高突变基因数，如果无txt，这部分的数量直接影响下游找DTG的数量，当然越多计算时间越长
- mem_gene2role 运算需要的内存，cpu已经内置设置为8


# Output
gene2role目录
  - 展示DTGs的距离分布图
  - 取距离最远的前5展示拓扑结构
  - 其它图和绘图数据


# Detail
> gene2role /opt/software/miniconda3/envs/gene2role/bin/
> louvain--01; gene2role--04
```shell
source /software/miniconda/bin/activate
conda info --envs
conda config --remove-key channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n gene2role r-base=4.3 python=3.12 -y
conda activate gene2role
pip install gensim #gensim 4.3.2 would require python >=3.12,<3.13.0a0 , which can be installed;
conda install conda-forge::r-seurat -y
conda install conda-forge::r-ggraph -y
conda install conda-forge::r-tidygraph -y
conda install conda-forge::r-tidyverse -y
conda install conda-forge::r-ggsignif -y
conda install anaconda::ipykernel -y
conda install conda-forge::r-optparse -y
conda install conda-forge::r-ggvenndiagram -y
conda install conda-forge::r-rpresto -y
conda install conda-forge::r-devtools -y

pip install futures
pip install fastdtw
pip install pandas
pip install matplotlib
pip install networkx
pip install python-louvain
```

# Reference
- [单细胞GRN新算法——识别『表达不变但功能巨变』的『差异拓扑』基因？首个角色驱动的基因嵌入算法~](https://mp.weixin.qq.com/s/3kGqRIhSuEWqS3Ykwhq2TA)
