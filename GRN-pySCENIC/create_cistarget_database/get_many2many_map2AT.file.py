# Title: get_many2many_map2AT.file.py
# Date: 20250524
# Coder: ydgenomics
# Description: Get one-to-one mapping and tf binding motif info
# Input: 
# Output: TF_motifs_many2many.tbl TF_gene_maped.txt
# Image: GRN-SCENIC-database--01 python
# Reference: kimi

import pandas as pd
from io import StringIO
from itertools import product
import argparse
parser = argparse.ArgumentParser(description="Get one-to-one mapping and tf binding motif info.")
parser.add_argument('--blastp_txt', default="/data/work/0.peanut/GRN/input/arabidopsis_db/blastp_results.txt", help="Input BLASTP results file")
parser.add_argument('--tf_list_txt', default="/data/work/0.peanut/GRN/input/Ath_TF_list.txt", help="Input Arabidopsis TF list file")
parser.add_argument('--tf_motifs_tbl', default="/data/work/0.peanut/GRN/output/Ath_TF_binding_motifs_information.tbl", help="Input TF binding motifs file")
parser.add_argument('--lowest_percent_identity', type=float, default=80, help="Lowest percent identity for filtering BLASTP results")
args = parser.parse_args()
input_blastp_txt = args.blastp_txt
input_AT_TF_list_txt = args.tf_list_txt
input_Ath_TF_motifs_tbl = args.tf_motifs_tbl
lowest_percent_identity = args.lowest_percent_identity


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

# Athough I will keep the many-to-many mapping, but I filter the low-quality mapping before.
# # # 过滤掉 percent_identity 小于 80 的行
blastp_df = blastp_df[blastp_df['percent_identity'] >= lowest_percent_identity] 

# # 按 query_id 分组，并在每个分组中选择 percent_identity 最高的行
# filtered_blastp_df = blastp_df.loc[blastp_df.groupby('query_id')['percent_identity'].idxmax()]
# print(filtered_blastp_df.head())
# len(filtered_blastp_df["query_id"])
# len(filtered_blastp_df["query_id"].unique())
# len(filtered_blastp_df["subject_id"].unique())

# # 按 subject_id 分组，并在每个分组中选择 percent_identity 最高的行
# filtered_blastp_df2 = filtered_blastp_df.loc[filtered_blastp_df.groupby('subject_id')['percent_identity'].idxmax()]
# print(filtered_blastp_df2.head())
# len(filtered_blastp_df2["query_id"])
# len(filtered_blastp_df2["query_id"].unique())
# len(filtered_blastp_df2["subject_id"].unique())

########## Get the AT name in subject_id2 ##########
# 读取TF_ID和Gene_ID的对应关系
with open(input_AT_TF_list_txt, 'r') as file:
    content = file.read()

tflist_df = pd.read_csv(StringIO(content), sep='\t', header=0)
tflist_df.head()
# 将 tflist_df 的 TF_ID 和 Gene_ID 转换为字典
tf_gene_dict = tflist_df.set_index('TF_ID')['Gene_ID'].to_dict()
# 在 filtered_blastp_df2 中匹配 subject_id 到 TF_ID，并添加 Gene_ID 列
blastp_df['subject_id2'] = blastp_df['subject_id'].map(tf_gene_dict)
blastp_df.head
# 去除没有匹配到的行
blastp_df = blastp_df.dropna(subset=['subject_id2'])
print(blastp_df.head())
len(blastp_df["query_id"])
len(blastp_df["query_id"].unique())
len(blastp_df["subject_id"].unique())
blastp_df[['query_id', 'subject_id', 'subject_id2']].head()

######### Get the TF binding motifs information ##########
data = pd.read_csv(input_Ath_TF_motifs_tbl, sep="\t")
#data = pd.read_csv("/data/work/0.peanut/GRN/input/ITAG4.1_MOTIF_PlantTFDB.tbl", sep="\t")
data.head()
len(data["motif_description"])
len(data["motif_description"].unique())
len(data["motif_name"])
len(data["motif_name"].unique())

# 创建一个空的 DataFrame 来存储最终结果
final_result = pd.DataFrame()
# 遍历 blastp_df 中的每个唯一 subject_id2
for subset in blastp_df['subject_id2'].unique():
    # 获取当前子集
    subset_df = blastp_df[blastp_df['subject_id2'] == subset]
    subset_data = data[data['motif_description'] == subset]
    # 如果 subset_df 或 subset_data 为空，跳过当前循环
    if subset_df.empty or subset_data.empty:
        continue
    # 创建一个空的列表来存储所有新行
    new_rows = []
    # 使用 itertools.product 生成笛卡尔积
    for (index, row), (_, data_row) in product(subset_df.iterrows(), subset_data.iterrows()):
        # 创建一个新的行，将 query_id 替换到 motif_description
        new_row = data_row.copy()
        new_row['motif_description'] = row['query_id']
        new_rows.append(new_row)
    # 将所有新行合并为一个 DataFrame
    result = pd.DataFrame(new_rows)
    # 将当前子集的结果追加到最终结果中
    final_result = pd.concat([final_result, result], ignore_index=True)

final_result['gene_name']=final_result['motif_description']
final_result[['motif_name', 'motif_description', 'gene_name']].head()
print(final_result)

final_result.to_csv("TF_motifs_many2many.tbl", sep="\t", index=False)

######### Get the TF gene information ##########
TF_list = final_result['motif_description'].unique()
# 将 NumPy 数组转换为 Pandas Series
TF_series = pd.Series(TF_list)
TF_series.to_csv("TF_gene_maped.txt", index=False, header=False, sep="\n")