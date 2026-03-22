import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

# ---
input_csv = 'Summary_Novelty.csv'
output_prefix = 'terpene'

df = pd.read_csv(input_csv)

cluster_counts = df['Louvain cluster'].value_counts()
valid_clusters = cluster_counts[cluster_counts >= 5].index
cluster_mapping = {old_name: f"Cluster {i+1}" for i, old_name in enumerate(valid_clusters)}
df_filtered = df[df['Louvain cluster'].isin(valid_clusters)].copy()
df_filtered['Cluster Name'] = df_filtered['Louvain cluster'].map(cluster_mapping)

def extract_number(cluster_name):
    return int(re.search(r'\d+', cluster_name).group())

name_order = sorted(df_filtered['Cluster Name'].unique(), key=extract_number)

#
np.random.seed(42)
mask = df_filtered['distance to closest MIBiG'] == 0.51
df_filtered.loc[mask, 'distance to closest MIBiG'] = np.random.normal(loc=0.51, scale=0.004, size=mask.sum())

#
custom_colors_hex = [
    "#008dff", "#ff7366", "#c701ff", "#4ecb8d", "#ff9d3a", 
    "#f9e858", "#d83034", "#8fd7d7", "#98c127", "#990099", "#004C99", "#994C00", "#FF3333"
]
color_map = {f"Cluster {i+1}": custom_colors_hex[i % len(custom_colors_hex)] for i in range(len(valid_clusters))}

# ==========================================
#
# ==========================================
sns.set_style("white")
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.edgecolor'] = 'black'

fig, ax = plt.subplots(figsize=(11, 6))

# 1.
sns.stripplot(
    data=df_filtered, 
    x='Cluster Name', 
    y='distance to closest MIBiG',
    order=name_order,
    palette=color_map,
    jitter=0.25,
    size=4.5,
    alpha=0.9,
    zorder=1,
    ax=ax,
    edgecolor='white',
    linewidth=0.4
)

# 2.
sns.boxplot(
    data=df_filtered, 
    x='Cluster Name', 
    y='distance to closest MIBiG',
    order=name_order,
    color="whitesmoke", 
    showfliers=False, 
    whis=(0, 100),  # <---
    width=0.5, 
    zorder=1,
    ax=ax,
    boxprops=dict(edgecolor='black', alpha=0.5, linewidth=0.5),
    medianprops=dict(color='black', linewidth=1.0),
    whiskerprops=dict(color='black', linewidth=1.0), 
    capprops=dict(color='black', linewidth=1.0)      
)

# 0.50 Cutoff
ax.axhline(y=0.50, color='red', linestyle='--', linewidth=0.5, zorder=0)

ax.set_ylim(-0.01, 0.55)
ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55])
ax.set_xlabel("Louvain Community", fontsize=12, fontweight='bold', labelpad=10)
ax.set_ylabel("Distance to closest MIBiG BGC\n(0.51 = Putatively Novel)", fontsize=12, fontweight='bold', labelpad=10)
plt.xticks(rotation=45, ha='right', fontsize=11, fontweight='bold')
plt.yticks(fontsize=11)
sns.despine(top=True, right=True)

plt.tight_layout()
plt.savefig(f"Panel_B_{output_prefix}_BoxPlot_FullWhiskers.png", dpi=300, bbox_inches='tight')
