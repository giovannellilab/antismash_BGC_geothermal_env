import pandas as pd
import numpy as np

#  ---
node_file = 'terpene_c007_nodetable.csv'
edge_file = 'terpene_c007_edgetable.csv'
gephi_prefix = 'terpene'
# ---------------------------------------

df_nodes = pd.read_csv(node_file)
df_edges = pd.read_csv(edge_file)

# 1. Edge Table
df_edges[['Source', 'Target']] = df_edges['name'].str.split(r' \(interacts with\) ', expand=True, regex=True)
gephi_edges = pd.DataFrame({'Source': df_edges['Source'], 'Target': df_edges['Target'], 'Type': 'Undirected'})
gephi_edges['Weight'] = 1.0 - df_edges['distance']
gephi_edges.to_csv(f'Gephi_Edges_{gephi_prefix}.csv', index=False)

# 2.
connected_nodes = set(gephi_edges['Source']).union(set(gephi_edges['Target']))
singletons = set(df_nodes['name']) - connected_nodes

# 3.
#
CUSTOM_PALETTE = [
    "#008dff", "#ff7366", "#c701ff", "#4ecb8d", "#ff9d3a", 
    "#f9e858", "#d83034", "#8fd7d7", "#98c127", "#990099", "#004C99", "#994C00", "#FF3333"
]

#
samples_only = df_nodes[~df_nodes['name'].str.startswith('BGC')]
cluster_sizes = samples_only['Louvain_Cluster'].value_counts()

cluster_map = {}
color_map = {}
for i, old_cluster in enumerate(cluster_sizes.index):
    new_name = f"Cluster {i+1}"
    cluster_map[old_cluster] = new_name
    color_map[new_name] = CUSTOM_PALETTE[i % len(CUSTOM_PALETTE)] # Palet biterse başa döner

def get_gephi_color(row):
    if str(row['name']).startswith('BGC'):
        return '#000000' # MIBiG'ler her zaman siyah
    new_c_name = cluster_map.get(row['Louvain_Cluster'], 'Other')
    return color_map.get(new_c_name, '#cccccc')

def get_new_cluster_name(row):
    if str(row['name']).startswith('BGC'):
        return 'MIBiG'
    return cluster_map.get(row['Louvain_Cluster'], 'Other')

# 4. Node Table
gephi_nodes = pd.DataFrame()
gephi_nodes['Id'] = df_nodes['name']
gephi_nodes['Label'] = df_nodes['name']
gephi_nodes['Cluster_Name'] = df_nodes.apply(get_new_cluster_name, axis=1)
gephi_nodes['Color'] = df_nodes.apply(get_gephi_color, axis=1)

# Singleton file creation
gephi_nodes = gephi_nodes[~gephi_nodes['Id'].isin(singletons)]
gephi_nodes.to_csv(f'Gephi_Nodes_{gephi_prefix}.csv', index=False)

