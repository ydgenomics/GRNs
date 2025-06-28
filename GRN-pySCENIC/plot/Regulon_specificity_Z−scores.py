# Title: Regulon_specificity_Z−scores.py
# Date: 20250525
# Coder: ydgenomics
# Description: plot the results of pySCENIC
# Input: 
# Output:
# Image: GRN-SCENIC-pySCENIC-database--01 python
# Reference: https://mp.weixin.qq.com/s/fTCrQdMuTy0eYNV6yh-Nwg

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import adjustText
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from adjustText import adjust_text
import argparse

parser = argparse.ArgumentParser(description="Plot regulon specificity Z-scores.")
parser.add_argument('--input_h5ad', type=str, default="/data/work/0.peanut/annotation/three_layers/H2014_dataget_Anno_rename_threelayers.h5ad", help='Path to input h5ad file')
parser.add_argument('--aucell_csv', type=str, default="/data/work/0.peanut/GRN/peanut/H2014/aucell.csv", help='Path to AUCell CSV file')
parser.add_argument('--cluster_key', type=str, default="cell", help='Key for cell type clustering in AnnData obs')
args = parser.parse_args()
input_h5ad = args.input_h5ad
aucell_csv = args.aucell_csv
cluster_key = args.cluster_key

# 读取数据
data2 = sc.read_h5ad(input_h5ad)
auc_mtx = pd.read_csv(aucell_csv, index_col = 'Cell')

rss_cellType = regulon_specificity_scores(auc_mtx, data2.obs[cluster_key])
#cats = sorted(list(set(data2.obs['cell'])))
cats = sorted(data2.obs[cluster_key].unique())

mat = rss_cellType.to_numpy()
row_mean = np.mean(mat, axis=1, keepdims=True)
row_std = np.std(mat, axis=1, keepdims=True)
norm_mat = (mat - row_mean) / row_std

rss_cellType_norm = pd.DataFrame(norm_mat, columns=rss_cellType.columns, index=rss_cellType.index)

# 确定子图的行数和列数
num_cats = len(cats)
cols = 4  # 每行的子图数量
rows = (num_cats + cols - 1) // cols  # 计算需要的行数

# 创建一个新的图形
fig = plt.figure(figsize=(20, rows * 5))  # 可以根据需要调整图形大小

# 遍历每个细胞类型和子图编号
for c, num in zip(cats, range(1, num_cats + 1)):
    x = rss_cellType_norm.T[c]
    ax = fig.add_subplot(rows, cols, num)  # 创建子图
    plot_rss(rss_cellType_norm, c, top_n=6, max_n=None, ax=ax)  # 绘制 RSS 数据
    ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05, x.max() + (x.max() - x.min()) * 0.05)  # 设置 y 轴范围
    for t in ax.texts:
        t.set_fontsize(8)  # 设置文本大小
    
    ax.set_ylabel('')  # 移除 y 轴标签
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-', color='black'), precision=0.001)  # 调整文本

# 隐藏多余的子图位置
for i in range(num_cats + 1, rows * cols + 1):
    fig.add_subplot(rows, cols, i)

fig.text(0.00, 0.5, 'Regulon specificity Z−scores in cell type', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()  # 自动调整子图参数，使之填充整个图像区域
plt.rcParams.update({
    'figure.autolayout': True,
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium'
})
pdf_file_path = 'Regulon_specificity_Z−scores.pdf'  # 指定保存的 PDF 文件路径
fig.savefig(pdf_file_path, bbox_inches='tight', pad_inches=0.1)  # 保存为 PDF
# 关闭图形，释放内存
plt.close(fig)