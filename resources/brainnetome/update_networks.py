import pandas as pd

# File paths
tsv_file = "brainnetome_network_mapping.csv"
text_file = "fsaverage.BN_Atlas.32k_fs_LR_updated.dlabel.labels.txt"
output_file = "fsaverage.BN_Atlas.32k_fs_LR_updated.dlabel.labels.withnetworks.txt"

# Load the TSV file into a DataFrame
mapping_df = pd.read_csv(tsv_file, sep="\,")

# Create a dictionary for quick lookup of New_Mapping
mapping_dict = dict(zip(mapping_df['Name'], mapping_df['New_Mapping']))



# Process the input text file and update labels
with open(text_file, 'r') as infile, open(output_file, 'w') as outfile:
    lines = infile.readlines()
    for i in range(0, len(lines), 2):  # Process lines in pairs (label and values)
        label = lines[i].strip()  # Get the label (e.g., L_A8m)
        values = lines[i + 1].strip()  # Get the values (e.g., 1 0 255 0 255)
        
        # Append the New_Mapping if the label exists in the mapping dictionary
        if label in mapping_dict:
            new_label = f"{label.replace('/', '-')}_{mapping_dict[label]}"  # Replace slashes with hyphens
        else:
            new_label = label  # Keep original label if no mapping found
        
        # Write updated label and values to the output file
        outfile.write(new_label + "\n")
        outfile.write(values + "\n")

print(f"Updated file saved to {output_file}")

