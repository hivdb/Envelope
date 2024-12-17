import csv
import re
import sys
import pandas as pd
import pickle

from hiv_seq_utils import read_fasta

pd.set_option('display.max_rows', 50)

# Constants
GENE = "GP41"

# Input files
hxb2_seq_path = "RefData/HXB2_" + GENE + "_AA.fasta"
mutations_file = "Mutations_" + GENE + ".xlsx"
(header, hxb2_seq) = read_fasta(hxb2_seq_path)
ref_seq_len = len(hxb2_seq)

# Output files
frequencies_all_output_file = f'{GENE}_Frequencies_All.xlsx'
frequencies_ge1pcnt_file = f'{GENE}_Frequencies_GE1pcnt.xlsx'
frequencies_horizontal_profile_file = f'{GENE}_Profile.xlsx'

# Read from the excel file which was created by gp120_align.py
mutations_df = pd.read_excel(mutations_file, header=[0], skiprows=[1])
print("mutations_df:\n", mutations_df)

subtype_pcnt = mutations_df['Subtype'].value_counts(normalize=True) * 100
#print(subtype_pcnt.head(10))
mutations_df = mutations_df.iloc[:, 3:]

def modify_cell(cell, pos):
    hxb2_ref = hxb2_seq[pos-1]
    cell = str(cell)
    if cell == '-':
        return hxb2_ref
    elif len(cell) >1:
        return 'INS'
    else: 
        return cell

for i, col in enumerate(mutations_df.columns, start=1):
    mutations_df[col] = mutations_df[col].apply(lambda cell: modify_cell(cell, i))
#print("Amino Acid Profile - Updated:\n", mutations_df, "\n")


frequencies = {}
for column in mutations_df.columns:
    value_counts = mutations_df[column].value_counts(normalize=True) * 100
    frequencies[column] = value_counts.round(1)
frequencies_df = pd.DataFrame(frequencies).fillna(0)
print("Frequencies dataframe\n", frequencies_df, "\n")


# Sort each column of frequencies_df in descending order of AA frequency
# Create two df, one with all frequencies and one with frequencies >= 1%
sorted_data = []
sorted_data_ge1pcnt = []
for column in frequencies_df.columns:
    # Get the sorted value frequencies (excluding zeros)
    sorted_values = frequencies_df[column][frequencies_df[column] > 0].sort_values(ascending=False)
    for value, frequency in sorted_values.items():
        sorted_data.append({
            "Pos": column,
            "AA": value,
            "Pcnt": frequency
        })
        if frequency < 1.0:
            continue
        sorted_data_ge1pcnt.append({
            "Pos": column,
            "AA": value,
            "Pcnt": frequency
        })

frequencies_sorted_df = pd.DataFrame(sorted_data)
frequencies_sorted_df.to_excel(frequencies_all_output_file, index=False)
frequencies_sorted_ge1pcnt_df = pd.DataFrame(sorted_data_ge1pcnt)
frequencies_sorted_ge1pcnt_df.to_excel(frequencies_ge1pcnt_file, index=False)

print("Parsed_frequencies >= 1%:\n", frequencies_sorted_df)

horiz_profile = {}
for mutations in sorted_data_ge1pcnt:
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
df.to_excel(frequencies_horizontal_profile_file, index=False)
#sys.exit("Exiting the program.")


