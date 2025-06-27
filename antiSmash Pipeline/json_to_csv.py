import json
import csv
import re
import os
import glob

# Define the base folder where the subfolders are located
base_folder = os.path.join(os.getcwd(), 'Pipeline Output')
# Find all JSON files in the subfolders of the base folder
json_files = glob.glob(os.path.join(base_folder, '*', '*.json'))

# Process each JSON file
for input_json_file in json_files:

    # Get the base name (without extension) to use for other files
    file_name = os.path.splitext(input_json_file)[0]

    # Construct output filenames using the base name
    output_csv_file = f"{file_name}.csv"

    # Load the JSON data
    with open(input_json_file, 'r') as json_file:
        data = json.load(json_file)

    data_list = []
    # Open CSV file for writing
    with open(output_csv_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the header row
        csv_writer.writerow(['Node ID', 'Product', 'Category'])

        # Extract and write relevant information from JSON to a CSV file
        for cluster in data['records']:
            if cluster['areas']:
                #Choose depending on the data from the json file
                # node_id = cluster['id']
                # node_id = re.search(r'NODE_(\d+)', node_id).group(1)
                # node_id = f"NODE {node_id}"
                node_id = cluster['id']
                node_id = re.search(r'contig_(\d+)', node_id).group(1)
                node_id = f"CONTIG {node_id}"
                # If there is more than one possible product for a node, choose the biggest one
                if len(cluster['areas']) > 1:
                    area_id = 0
                    idx_dict = None
                    for index, x in enumerate(cluster['areas']):
                        value = x['end'] - x['start']
                        if area_id < value:
                            area_id = value
                            idx_dict = index
                    product = cluster['areas'][idx_dict]['protoclusters']['0']['product']
                    category = cluster['areas'][idx_dict]['protoclusters']['0']['category']
                else:
                    product = cluster['areas'][0]['protoclusters']['0']['product']
                    category = cluster['areas'][0]['protoclusters']['0']['category']
                    data = cluster['seq']

                # Write the row to CSV
                csv_writer.writerow([node_id, product, category])

                # Create a dictionary for the current cluster
                cluster_info = {
                    'node_id': node_id,
                    'product': product,
                    'category': category
                }

                # Append the dictionary to the data list
                data_list.append(cluster_info)

    print(f'CSV file has been created: {output_csv_file}')

