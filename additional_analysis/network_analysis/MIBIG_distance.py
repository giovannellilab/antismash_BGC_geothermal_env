# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 15:14:34 2026


"""

import pandas as pd
import numpy as np


# 1.
node_file = 'terpene_c007_nodetable.csv'
edge_file = 'terpene_c007_edgetable.csv' #

df_nodes = pd.read_csv(node_file)
df_edges = pd.read_csv(edge_file)

#
df_edges[['Source', 'Target']] = df_edges['name'].str.split(r' \(interacts with\) ', expand=True, regex=True)

#
bgc_edges = df_edges[df_edges['Source'].str.startswith('BGC') | df_edges['Target'].str.startswith('BGC')]

# 2.
records = []
for _, row in bgc_edges.iterrows():
    if row['Source'].startswith('BGC') and not row['Target'].startswith('BGC'):
        records.append({'BGC ID': row['Target'], 'closest MIBiG BGC': row['Source'], 'distance to closest MIBiG': float(row['distance'])})
    elif row['Target'].startswith('BGC') and not row['Source'].startswith('BGC'):
        records.append({'BGC ID': row['Source'], 'closest MIBiG BGC': row['Target'], 'distance to closest MIBiG': float(row['distance'])})

dist_df = pd.DataFrame(records)

#
if not dist_df.empty:
    #
    dist_df = dist_df.sort_values('distance to closest MIBiG').drop_duplicates('BGC ID', keep='first')
else:
    dist_df = pd.DataFrame(columns=['BGC ID', 'closest MIBiG BGC', 'distance to closest MIBiG'])

# 3.
df_samples = df_nodes[~df_nodes['name'].str.startswith('BGC')].copy()
df_samples = df_samples.rename(columns={'name': 'BGC ID', 'Louvain_Cluster': 'Louvain cluster'})

#
final_df = pd.merge(df_samples[['BGC ID', 'Louvain cluster']], dist_df, on='BGC ID', how='left')

# 4.
#
final_df['Is_0.51_Assigned'] = final_df['distance to closest MIBiG'].isna().map({True: 'Yes (No MIBiG connection)', False: 'No (Directly observed)'})
final_df['closest MIBiG BGC'] = final_df['closest MIBiG BGC'].fillna('None')
final_df['distance to closest MIBiG'] = final_df['distance to closest MIBiG'].fillna(0.51)

#
columns_order = ['BGC ID', 'Louvain cluster', 'closest MIBiG BGC', 'distance to closest MIBiG', 'Is_0.51_Assigned']
final_df = final_df[columns_order]

output_filename = 'Ana_Summary_Novelty_terpene.csv'
final_df.to_csv(output_filename, index=False)
