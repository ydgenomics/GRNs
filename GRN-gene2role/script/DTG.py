### Date: 250729
### gene2role

import pandas as pd
import numpy as np
import os
import time
print(f'Now work path: {os.getcwd()}')
os.chdir('/Gene2Role/reproduce/')
print(f'Now work path: {os.getcwd()}')
import sys
sys.path.insert(0, '/Gene2Role/reproduce')   # 放在 import dtg_utils 之前
from dtg_utils import *
from draw_1_hop import *
import argparse

parser = argparse.ArgumentParser(description='Process glioblastoma data.')
parser.add_argument('--input_emb', type=str, default='/data/work/Single-Cell-Pipeline/output/GRN-gene2role/all_cell/result/sample.emb', help='Path to embedding file')
parser.add_argument('--index_tsv', type=str, default='/data/work/Single-Cell-Pipeline/output/GRN-gene2role/all_cell/result/splitMatrix/index_tracker.tsv', help='Path to index TSV file')
parser.add_argument('--input_edgelist', type=str, default='/data/work/Single-Cell-Pipeline/output/GRN-gene2role/all_cell/result/sample_merged.edgelist', help='Path to edgelist file')
parser.add_argument('--work_path', type=str, default='/data/work/Single-Cell-Pipeline/output/GRN-gene2role/all_cell/result/tools/SignedS2V/graph', help='Path to edgelist file')
parser.add_argument('--n_top', type=int, default=5, help='Number of top nodes to plot')
args = parser.parse_args()
input_emb = args.input_emb
index_tsv = args.index_tsv
input_edgelist = args.input_edgelist
work_path = args.work_path
n_top = args.n_top

os.chdir(work_path)
glioblastoma_embedding = pd.read_csv(input_emb,sep=' ',skiprows=1,header=None,index_col=0)
glioblastoma_index = pd.read_csv(index_tsv,sep='\t')
if 'Unnamed: 0' not in glioblastoma_index.columns:
    glioblastoma_index.insert(
        loc=0,                       # 插到第一列
        column='Unnamed: 0',
        value=glioblastoma_index.index
    )
    # 1. 去掉无意义的重复行名（让索引变为 0,1,2...）
    glioblastoma_index.reset_index(drop=True, inplace=True)
    # 2. 把除第一列外的所有列转成 int64
    cols_to_int = glioblastoma_index.columns.difference(['Unnamed: 0'])
    glioblastoma_index[cols_to_int] = (
        glioblastoma_index[cols_to_int]
        .apply(pd.to_numeric, errors='coerce')   # 如有异常值先转成 NaN
        .fillna(0)                              # 可选：NaN 用 0 补
        .astype('int64')
    )
    # 3. 检查结果
    print(glioblastoma_index.head())
    print(glioblastoma_index.dtypes)   # 确认 dtype

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

# time.sleep(10)   # 暂停 10 秒

glioblastoma_distance.to_csv('glb_eeisp_distance.csv', sep='\t', index=False)

### plot the nodes with largest distance
# read the orginal edgelist
glioblastoma_edgelist = input_edgelist
edgelist = pd.read_csv(glioblastoma_edgelist, sep='\t', header=None)
print(edgelist.head())
edgelist_list = edgelist.iloc[:, :2].stack().unique().tolist()
print(edgelist_list[0:5])


if 'distance_avg' in glioblastoma_distance.columns:
    plot_avg_sd(glioblastoma_distance, out_dir=".")
    ### multi
    # ---- avg -----
    output_path = os.getcwd()+'/avg/'
    os.makedirs("avg", exist_ok=True)
    MEP_distance = glioblastoma_distance
    MEP_distance_sort_mean = MEP_distance.sort_values('distance_avg', ascending = False)
    n_cluster = len(glioblastoma_distance.columns) - 2
    rows_index_list = [row[0:n_cluster].tolist() for index, row in MEP_distance_sort_mean.head(n_top).iterrows()]
    for index_list in rows_index_list:
        # index_list = [int(x) for x in index_list]
        # index_list = [str(x) for x in index_list]
        # edgelist_list = [str(x) for x in edgelist_list]
        batch_plot_1_hop(glioblastoma_edgelist,index_list[1:],output_path,index_list[0], edgelist_list)
    # ---- sd -----
    output_path = os.getcwd()+'/sd/'
    os.makedirs("sd", exist_ok=True)
    MEP_distance_sort_sd = MEP_distance.sort_values('distance_sd', ascending = False)
    n_cluster = len(glioblastoma_distance.columns) - 2
    rows_index_list = [row[0:n_cluster].tolist() for index, row in MEP_distance_sort_sd.head(n_top).iterrows()]
    for index_list in rows_index_list:
        # index_list = [int(x) for x in index_list]
        # index_list = [str(x) for x in index_list]
        batch_plot_1_hop(glioblastoma_edgelist,index_list[1:],output_path,index_list[0], edgelist_list)
else: 
    plot_frequencies(glioblastoma_distance, out_dir=".")
    ### two group
    path = os.getcwd()+'/'
    id_list = glioblastoma_distance.iloc[:n_top, [0, 1, 2]].values.tolist()
    print(id_list)
    for index_list in id_list:
        print(index_list[0])
        batch_plot_1_hop(glioblastoma_edgelist,index_list[1:],path,index_list[0], edgelist_list)


### plot the node which is our target
# glioblastoma_distance = pd.read_csv('/data/work/Single-Cell-Pipeline/output/GRN-gene2role/plot/glb_eeisp_distance.csv', sep='\t')
# egr = glioblastoma_distance[glioblastoma_distance['Unnamed: 0'] == 'Ga07g01598'].values.tolist()[0]
# egr
# batch_plot_1_hop(glioblastoma_edgelist,egr[1:],path,egr[0], edgelist_list)