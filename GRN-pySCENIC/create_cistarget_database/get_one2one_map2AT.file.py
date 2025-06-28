# Title: get_one2one_map2AT.file.py
# Date: 20250523
# Coder: ydgenomics
# Description: Get one-to-one mapping and tf binding motif info
# Input: 
# Output: TF_motifs_maped.tbl TF_gene_maped.txt
# Image: GRN-SCENIC-database--01 python
# Reference: kimi

import pandas as pd
from io import StringIO
import argparse
parser = argparse.ArgumentParser(description="Get one-to-one mapping and tf binding motif info.")
parser.add_argument('--blastp_txt', default="/data/work/0.peanut/GRN/input/arabidopsis_db/blastp_results.txt", help="Input BLASTP results file")
parser.add_argument('--tf_list_txt', default="/data/work/0.peanut/GRN/input/Ath_TF_list.txt", help="Input Arabidopsis TF list file")
parser.add_argument('--tf_motifs_tbl', default="/data/work/0.peanut/GRN/output/Ath_TF_binding_motifs_information.tbl", help="Input TF binding motifs file")
args = parser.parse_args()
input_blastp_txt = args.blastp_txt
input_AT_TF_list_txt = args.tf_list_txt
input_Ath_TF_motifs_tbl = args.tf_motifs_tbl

########## Get the one-to-one mapping result for blastp ##########
with open(input_blastp_txt, 'r') as file:
    content = file.read()

# 将内容转换为DataFrame
blastp_df = pd.read_csv(StringIO(content), sep='\t', header=None)
# 添加列名
blastp_df.columns = [
    'query_id', 'subject_id', 'percent_identity', 'alignment_length',
    'mismatches', 'gap_opens', 'query_start', 'query_end',
    'subject_start', 'subject_end', 'e_value', 'bit_score'
]
# 查看前几行
print(blastp_df.head())
len(blastp_df["query_id"])
len(blastp_df["query_id"].unique())
len(blastp_df["subject_id"].unique())

# 按 query_id 分组，并在每个分组中选择 percent_identity 最高的行
filtered_blastp_df = blastp_df.loc[blastp_df.groupby('query_id')['percent_identity'].idxmax()]
print(filtered_blastp_df.head())
len(filtered_blastp_df["query_id"])
len(filtered_blastp_df["query_id"].unique())
len(filtered_blastp_df["subject_id"].unique())

# 按 subject_id 分组，并在每个分组中选择 percent_identity 最高的行
filtered_blastp_df2 = filtered_blastp_df.loc[filtered_blastp_df.groupby('subject_id')['percent_identity'].idxmax()]
print(filtered_blastp_df2.head())
len(filtered_blastp_df2["query_id"])
len(filtered_blastp_df2["query_id"].unique())
len(filtered_blastp_df2["subject_id"].unique())

########## Get the AT name in subject_id2 ##########
# 读取TF_ID和Gene_ID的对应关系
with open(input_AT_TF_list_txt, 'r') as file:
    content = file.read()

tflist_df = pd.read_csv(StringIO(content), sep='\t', header=0)
tflist_df.head()
# 将 tflist_df 的 TF_ID 和 Gene_ID 转换为字典
tf_gene_dict = tflist_df.set_index('TF_ID')['Gene_ID'].to_dict()
# 在 filtered_blastp_df2 中匹配 subject_id 到 TF_ID，并添加 Gene_ID 列
filtered_blastp_df2['subject_id2'] = filtered_blastp_df2['subject_id'].map(tf_gene_dict)
# 去除没有匹配到的行
filtered_blastp_df2 = filtered_blastp_df2.dropna(subset=['subject_id2'])
print(filtered_blastp_df2.head())
len(filtered_blastp_df2["query_id"])
len(filtered_blastp_df2["query_id"].unique())
len(filtered_blastp_df2["subject_id"].unique())

######### Get the TF binding motifs information ##########
data = pd.read_csv(input_Ath_TF_motifs_tbl, sep="\t")
# 创建一个字典，将 filtered_blastp_df2 中的 subject_id2 映射到 query_id
blastp_dict = dict(zip(filtered_blastp_df2['subject_id2'], filtered_blastp_df2['query_id']))
# 在 data 中匹配 subject_id2 到 query_id
data['motif_description'] = data['motif_description'].map(blastp_dict)
data['gene_name'] = data['gene_name'].map(blastp_dict)
# 删除没有匹配到的行
data = data.dropna(subset=['motif_description', 'gene_name'])
# 查看结果
print(data.head())
data.to_csv("TF_motifs_maped.tbl", sep="\t", index=False)

######### Get the TF gene information ##########
query_ids = filtered_blastp_df2['query_id']
# 保存为文本文件，每个值占一行
output_file = "TF_gene_maped.txt"
query_ids.to_csv(output_file, index=False, header=False, sep="\n")