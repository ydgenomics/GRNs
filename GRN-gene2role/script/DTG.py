import pandas as pd
import numpy as np
import os
os.chdir('/Gene2Role/reproduce/')
from dtg_utils import *
from draw_1_hop import *
import argparse

glioblastoma_embedding = pd.read_csv('/data/work/GRN-gene2role/plot/rep_data/rep_data/result2/human_glioblastoma/0_12_hours.emb',
                                    sep=' ',
                                    skiprows=1,
                                    header=None,
                                    index_col=0)
glioblastoma_index = pd.read_csv('/data/work/GRN-gene2role/test3/result/splitMatrix/index_tracker.tsv',
                                sep='\t')
glioblastoma_index.head(5)

glioblastoma_distance = cal_node_distances(glioblastoma_index, glioblastoma_embedding)

def plot_avg_sd(df,
                avg='distance_avg',
                sd='distance_sd',
                out_dir='.',
                out_name='avg_sd_plot.pdf'):
    """
    绘制 average-vs-SD 散点图并保存为 PDF。
    参数
    ----
    df      : DataFrame
    avg     : 平均距离列名
    sd      : 标准差列名
    out_dir : 输出目录（不存在会自动创建）
    out_name: 输出文件名，必须以 .pdf 结尾
    """
    os.makedirs(out_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 5))
    df.plot(kind='scatter', x=avg, y=sd,
            color='blue', alpha=0.5, ax=ax)
    ax.set_xlabel('Average Distance')
    ax.set_ylabel('SD of Distance')

    out_path = os.path.join(out_dir, out_name)
    fig.savefig(out_path, format='pdf', bbox_inches='tight')
    plt.close(fig)          # 避免后续 cell 里重复显示
    print(f"Saved: {out_path}")

def plot_frequencies(df,
                     column_names='distance',
                     out_dir='.',
                     out_name='distance_hist.pdf'):
    """
    绘制指定列的直方图并保存为 PDF。
    参数
    ----
    df        : DataFrame
    column_names: 要画直方图的列名（单个字符串）
    out_dir   : 输出目录（不存在会自动创建）
    out_name  : 输出文件名，必须以 .pdf 结尾
    """
    os.makedirs(out_dir, exist_ok=True)

    data = df[column_names].dropna()          # 去掉 NaN，避免空图报错
    if data.empty:
        print("Warning: no valid data to plot.")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(data, bins=30, edgecolor='k')
    ax.set_xlabel('Distance')
    ax.set_ylabel('Frequency')

    out_path = os.path.join(out_dir, out_name)
    fig.savefig(out_path, format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")

plot_frequencies(glioblastoma_distance, out_dir="/data/work/GRN-gene2role/test2/plot")

glioblastoma_distance.to_csv('/data/work/GRN-gene2role/test3/plot/glb_eeisp_distance.csv', sep='\t', index=False)

### plot the nodes with largest distance
# read the orginal edgelist
glioblastoma_edgelist = '/data/work/GRN-gene2role/plot/rep_data/rep_data/result2/human_glioblastoma/0_12_hours_merged.edgelist'
edgelist = pd.read_csv(glioblastoma_edgelist, sep='\t', header=None)
print(edgelist.head())
edgelist_list = edgelist.iloc[:, :2].stack().unique().tolist()
print(edgelist_list[0:5])
# save path
path = '/data/work/GRN-gene2role/test3/plot/'

id_list = glioblastoma_distance.iloc[:5, [0, 1, 2]].values.tolist()
print(id_list)
for index_list in id_list:
    print(index_list[0])
    batch_plot_1_hop(glioblastoma_edgelist,index_list[1:],path,index_list[0], edgelist_list)

### plot the node which is our target
egr = glioblastoma_distance[glioblastoma_distance['Unnamed: 0'] == 'GMNN'].values.tolist()[0]
egr
batch_plot_1_hop(glioblastoma_edgelist,egr[1:],path,egr[0], edgelist_list)