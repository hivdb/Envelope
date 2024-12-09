import csv
import re
import pandas as pd
import pickle

from hiv_seq_utils import read_fasta

pd.set_option('display.max_rows', 50)
hxb2 = "RefData/HXB2_GP120_AA.fasta"
(header, hxb2_seq) = read_fasta(hxb2)

# Once the dataframe has been serialized I commented out the following 2 
# blocks of code
# Read from the excel file which was created by gp120_align.py
mutations_file = "mutations.csv"
mutations_df = pd.read_csv(mutations_file, header=[0], skiprows=[1])
print(mutations_df)

subtype_pcnt = mutations_df['Subtype'].value_counts(normalize=True) * 100
print(subtype_pcnt.head(10))
subtypes = ['B', 'C', '01_AE', 'A1', '02_AG', 'D', 'G']

def modify_cell(cell, pos):
    hxb2_ref = hxb2_seq[pos-1]
    cell = str(cell)
    if cell == '-':
        return hxb2_ref
    elif len(cell) >1:
        return 'INS'
    else: 
        return cell

mutations_df = mutations_df.iloc[:, 3:]
print("Amino Acid Profile:\n", mutations_df)

for i, col in enumerate(mutations_df.columns, start=1):
    mutations_df[col] = mutations_df[col].apply(lambda cell: modify_cell(cell, i))
print("Amino Acid Profile - Updated:\n", mutations_df)

frequencies = {}
for column in mutations_df.columns:
    value_counts = mutations_df[column].value_counts(normalize=True) * 100
    frequencies[column] = value_counts.round(1)

with open('Frequencies.txt', 'w') as f:
   print(frequencies, file=f)

# Convert the frequencies dictionary to a DataFrame
frequencies_df = pd.DataFrame(frequencies).fillna(0)

# Step 2: Display the frequencies DataFrame
print("Frequencies dataframe\n", frequencies_df)

# Parse the DataFrame
parsed_data = []
for column in frequencies_df.columns:
    # Get the sorted value frequencies (excluding zeros)
    sorted_values = frequencies_df[column][frequencies_df[column] > 0].sort_values(ascending=False)
    for value, frequency in sorted_values.items():
        if frequency < 1.0:
            continue
        parsed_data.append({
            "Pos": column,
            "AA": value,
            "Pcnt": frequency
        })

frequencies_parsed_df = pd.DataFrame(parsed_data)
print("Parsed_frequencies >= 1%", frequencies_parsed_df)
frequencies_parsed_df.to_excel("Profile1.xlsx", index=False)

horiz_profile = {}
for mutations in parsed_data:
    position = int(mutations['Pos'][1:])
    aa = mutations['AA']
    pcnt = round(mutations['Pcnt'])
    if position not in horiz_profile:
        horiz_profile[position] = []
    horiz_profile[position].append({'AA': aa, 'Pcnt': pcnt})

# Sort and prepare the data for DataFrame
# Initialize empty lists for each position
data = {pos: [] for pos in sorted(horiz_profile.keys())}  

# Determine the maximum number of mutations per position
max_length = max(len(values) for values in horiz_profile.values())  

for pos in sorted(horiz_profile.keys()):
    # Sort by 'Pcnt' in descending order and format the entries
    sorted_entries = sorted(horiz_profile[pos], key=lambda x: x['Pcnt'], reverse=True)
    formatted_entries = ['{} {}%'.format(entry['AA'], entry['Pcnt']) for entry in sorted_entries]
    
    # Pad the list if it's shorter than max_length
    while len(formatted_entries) < max_length:
        formatted_entries.append('')
    
    # Assign the formatted list to the correct position
    data[pos] = formatted_entries

df = pd.DataFrame(data)
print(df)
df.to_excel('horiz_profile.xlsx', index=False)


