# Gene2Role
1. network construction [pearseam/eeisp]
2. embedding generation [Struc2vec]
3. downstream analysis [plot_and_intepretation]

```shell
python pipeline.py --help
usage: pipeline.py [-h] [--project PROJECT] [--reindex] [--output [OUTPUT]] [--dimensions DIMENSIONS] [--walk-length WALK_LENGTH] [--num-walks NUM_WALKS] [--window-size WINDOW_SIZE]
                   [--until-layer UNTIL_LAYER] [--iter ITER] [--workers WORKERS] [--OPT1] [--OPT2] [--OPT3] [--scalefree] [--correlation_threshold CORRELATION_THRESHOLD] [--matrix MATRIX]
                   [--mode MODE] [--threCDI THRECDI] [--threEEI THREEEI] [--cell_metadata CELL_METADATA]
                   TaskMode CellType EmbeddingMode input [input ...]

Gene2Role pipeline script

positional arguments:
  TaskMode              Task mode. 1: run SignedS2V for an edgelist file. 2: run spearman and SignedS2V from gene X cell count matrix. 3: run eeisp and SignedS2V from gene X cell count matrix.
  CellType              Cell type. 1: single cell-type. 2: multiple cell-type in one matrix. 3: multiple gene X cell count matrixes.
  EmbeddingMode         Embedding mode. 1: single network embedding. 2: multiple network embedding, only work if the previous argument is 2.
  input

options:
  -h, --help            show this help message and exit
  --project PROJECT     Project name which will be used as the folder name.
  --reindex             Flag for reindexing genes. Needed if gene name contained in input file.
  --output [OUTPUT]     [SignedS2V] Output emb path, if Not given, follow input file name
  --dimensions DIMENSIONS
                        [SignedS2V] Number of dimensions. Default is 128.
  --walk-length WALK_LENGTH
                        [SignedS2V] Length of walk per source. Default is 80.
  --num-walks NUM_WALKS
                        [SignedS2V] Number of walks per source. Default is 100.
  --window-size WINDOW_SIZE
                        [SignedS2V] Context size for optimization. Default is 10.
  --until-layer UNTIL_LAYER
                        [SignedS2V] Calculation until the layer. Default is 3.
  --iter ITER           [SignedS2V] Number of epochs in SGD
  --workers WORKERS     Number of parallel workers. Default is 8.
  --OPT1                optimization 1
  --OPT2                optimization 2
  --OPT3                optimization 3
  --scalefree           scale free flag
  --correlation_threshold CORRELATION_THRESHOLD
                        [spearman] Threshold for filtering high correlations. (default:0.4)
  --matrix MATRIX       [ESSIP/spearman/split cell] CSV file for gene-cell matrix.
  --mode MODE           [ESSIP] mode for (1) percentage threshold or (2) solid threshold for CDI and EEI (default: 1)
  --threCDI THRECDI     [ESSIP] Threshold for CDI (default: 0.5). When mode = 1, filter out top value*100%; when mode = 1, filter out those >= value.
  --threEEI THREEEI     [ESSIP] Threshold for EEI (default: 0.5). When mode = 1, filter out top value*100%; when mode = 1, filter out those >= value.
  --cell_metadata CELL_METADATA
                        [split cell] CSV file with cell type information.
```

[2025(BMC bioinformatics)_Gene2role: a role-based gene embedding method for comparative analysis of signed gene regulatory networks]()

# 环境配置 
> gene2role /opt/software/miniconda3/envs/gene2role/bin/
```shell
source /opt/software/miniconda3/bin/activate
conda create -n gene2role r-base=4.3 python=3.12 -y
conda activate gene2role
pip install gensim #gensim 4.3.2 would require python >=3.12,<3.13.0a0 , which can be installed;
conda install conda-forge::r-seurat -y
pip install futures
pip install fastdtw
pip install pandas
pip install matplotlib
```
```R
library(Seurat)
library(data.table)
```
```python
# conda install python=3.13
import subprocess
import argparse
import os
import pandas as pd
from logging import getLogger, INFO, DEBUG
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed # pip install futures
import multiprocessing
import copy
import gensim.models
```

[CellOracle](https://www.nature.com/articles/s41586-022-05688-9) integrates scATAC-seq and scRNA-seq data, leveraging transcription factor (TF) binding motifs and co-expression information to infer GRNs, providing a more detailed and comprehensive understanding of gene regulatory mechanisms.