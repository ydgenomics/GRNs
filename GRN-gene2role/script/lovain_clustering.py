# louvain--01
print("--------------------------- Lovain_clustering ---------------------------")
import networkx as nx
import community as community_louvain
import pandas as pd
import numpy as np
import os
os.chdir("/data/work/Single-Cell-Pipeline/output/GRN-gene2role/test")
# 谁作为anchor，就提交哪一个子edgelist
dir_list = {
    'test1':'/data/users/yangdong/yangdong_db936836d4034b638f5c86db02932db1/online/GRN-gene2role/plot/rep_data/rep_data/result2/human_glioblastoma/0_12_hours_merged.edgelist',
    'test2':'/data/users/yangdong/yangdong_db936836d4034b638f5c86db02932db1/online/GRN-gene2role/plot/rep_data/rep_data/result2/GMP_CD14_Monocytes/GMP_CD14_Monocytes_merged.edgelist',
}
louvain_res=0.8
for project, dir_ in dir_list.items():
    df = pd.read_csv(dir_, sep='\t', header=None)
    df = df[df[2] > 0]
    G = nx.from_pandas_edgelist(df,0,1,edge_attr=None, create_using=nx.Graph)
    connected_components = nx.connected_components(G)
    # Find the largest connected component
    largest_component = max(connected_components, key=len)
    # Create a new graph containing only the largest component
    largest_component_graph = G.subgraph(largest_component)
    connected_components = nx.connected_components(G)
    for x in connected_components:
        print (x)
    # Perform community detection using the Louvain method
    partition = community_louvain.best_partition(largest_component_graph, resolution=louvain_res)
    print (community_louvain.modularity(partition, largest_component_graph))
    df_out = pd.DataFrame(partition, index=['community']).T
    df_out.insert(0, "gene", df_out.index)
    print (df_out['community'].unique())
    df_out.to_csv(f"Louvain_{project}.csv", index=None, header=True, sep='\t')
    # Print the detected communities
    # print("Communities:")
    # for node, community in partition.items():
    #     print(f"Node {node} belongs to community {community}")

print("----------------------------- gene module analysis ----------------------------")
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/Gene2Role/reproduce')   # 放在 import dtg_utils 之前
from dtg_utils import *

def cal_pair_distance(row, embedding):
	# Extract indices
	index1, index2 = row[-2], row[-1]
	# Check if indices exist in embedding
	if index1 in embedding.index and index2 in embedding.index:
		return euclidean(embedding.loc[index1], embedding.loc[index2])
	else:
		return np.nan
    
def cal_batch_pair_distance(embedding, embedding_index, cluster_df, base_index_col = 1):
    cluster_df['distance_change'] = np.nan
    for index, row in cluster_df.iterrows():
        based_index = int(row['gene'])
        match_row = embedding_index[embedding_index.iloc[:, base_index_col] == based_index]
        match_row = match_row.values.tolist()[0]
        if match_row:
            cluster_df.at[index, 'distance_change'] = cal_pair_distance(match_row, embedding)
        else:
            print('Is there anything wrong with your input data?')
    return cluster_df

glb_embedding = pd.read_csv('/data/work/GRN-gene2role/plot/rep_data/rep_data/result2/human_glioblastoma/0_12_hours.emb',sep=' ',skiprows=1,header=None,index_col=0)
glb_index = pd.read_csv('/data/work/GRN-gene2role/test3/result/splitMatrix/index_tracker.tsv',sep='\t')
glb_index.head(5)

glb_lovain_0 = pd.read_csv('/data/work/Single-Cell-Pipeline/output/GRN-gene2role/test/Louvain_test1.csv',sep='\t')
glb_lovain_0['community'].value_counts()

glb_lovain_distance_0 = cal_batch_pair_distance(glb_embedding, glb_index, glb_lovain_0)
glb_distance_info_0 = glb_lovain_distance_0.groupby('community')['distance_change'].agg(['mean', 'std', 'count', 'size'])
glb_distance_info_0['nan_ratio'] = (glb_distance_info_0['size'] - glb_distance_info_0['count']) / glb_distance_info_0['size']
glb_distance_info_0.head(5)

glb_distance_info_0
glb_distance_info_0.to_csv('glb_distance_0.csv', sep='\t',quoting=False)