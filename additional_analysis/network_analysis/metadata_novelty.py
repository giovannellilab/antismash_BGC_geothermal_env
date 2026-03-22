import pandas as pd
import networkx as nx
import os

def get_mapping_info(node_id):
    """نقشه‌برداری گره‌ها به رنگ و شکل طبق دستورالعمل مشتری"""
    node_str = str(node_id)
    if node_str.startswith('BGC'):
        return "MIBiG (Known)", "#000000", "diamond"

    mapping = {
        "ARG": ("SA-CVZ (backarc)", "#E69F00"),
        "CH20": ("SA-CVZ (volcanic arc)", "#009E73"),
        "CHL22": ("SA-CVZ (volcanic arc)", "#009E73"),
        "CR": ("CAVA", "#56B4E9"),
        "ISL21": ("RVB", "#F0E442"),
        "MNG23": ("CAOB", "#CC79A7"),
        "TR_": ("CAOB", "#CC79A7"),
        "SJ_": ("CAOB", "#CC79A7"),
        "NR_": ("CAOB", "#CC79A7"),
        "KJ_": ("CAOB", "#CC79A7"),
        "HYD22": ("TLVP", "#D55E00"),
        "AEO": ("AAVP", "#0072B2"),
        "FEN": ("Campania", "#725CFF"),
        "GUY": ("Guaymas Basin", "#1AC9A3"),
        "EPR": ("EPR", "#017371")
    }
    
    for prefix, (label, color) in mapping.items():
        if node_str.startswith(prefix):
            return label, color, "ellipse"
    
    return "Other/Unknown", "#D3D3D3", "ellipse"

def process_bgc_data(file_dict):
    all_edges = []
    all_node_metadata = []
    summary_report = []

    for file_path, category in file_dict.items():
        if not os.path.exists(file_path):
            print(f"Skipping {category}: File {file_path} not found.")
            continue
        
        print(f"Processing {category}...")
        # 1.
        df = pd.read_csv(file_path, sep='\t')
        all_edges.append(df)
        
        # 2 NetworkX
        G = nx.Graph()
        for _, row in df.iterrows():
            G.add_edge(row['Record_a'], row['Record_b'])
        
        families = list(nx.connected_components(G))
        total_fams = len(families)
        known_fams = 0
        
        # 3.
        category_nodes = []
        for fam in families:
            is_known = any(str(node).startswith('BGC') for node in fam)
            if is_known:
                known_fams += 1
            
            for node in fam:
                label, color, shape = get_mapping_info(node)
                node_info = {
                    'NodeID': node,
                    'Category': category,
                    'Region': label,
                    'Hex_Color': color,
                    'Node_Shape': shape,
                    'Is_Known_Family': "Yes" if is_known else "No"
                }
                category_nodes.append(node_info)
                all_node_metadata.append(node_info)
        
        #
        pd.DataFrame(category_nodes).drop_duplicates(subset='NodeID').to_csv(f"{category}_metadata.csv", index=False)
        
        #
        summary_report.append({
            'Category': category,
            'Total_Families': total_fams,
            'Known_Families': known_fams,
            'Novel_Families': total_fams - known_fams,
            'Novelty_Percentage': f"{(1 - known_fams/total_fams)*100:.2f}%" if total_fams > 0 else "0%"
        })

    # 4 (Combined)
    if all_edges:
        pd.concat(all_edges).to_csv("Combined_Network.txt", sep='\t', index=False)
        pd.DataFrame(all_node_metadata).drop_duplicates(subset=['NodeID', 'Category']).to_csv("Combined_Metadata.csv", index=False)
    
    return pd.DataFrame(summary_report)

#
my_files = {
    'NRPS_c0.5.txt': 'NRPS',
    'PKS_c0.5.txt': 'PKS',
    'RiPP_c0.5.txt': 'RiPP',
    'terpene_c0.5.txt': 'Terpene'
}

#
final_report = process_bgc_data(my_files)
print("\n--- Novelty Analysis Report ---")
print(final_report)
final_report.to_csv("Final_Novelty_Report.csv", index=False)