import csv
import re
import sys
import os
import pandas as pd

from hiv_seq_utils import (read_fasta, 
                           align_2_aa_seqs, 
                           compare_2_strings, 
                           check_hxb2_start)

pd.set_option('display.max_rows', 50)
GENE = "GP120"

# The LANL_FilteredAASeqs file has amino acid sequences downloaded from GenBank
# that were in the original file LANL_HIV1_FLT_2021_accessions.txt
# The filtered file has only those sequences with more than 450 amino acids
env_sequences_file = "RefData/LANL_FilteredAASeqs_Dec15.xlsx"

# This file has the subtypes from LANL
#metadata_file = "RefData/LANL_gp120_search_1pp_13075_w_subtypes_June22.xlsx"
metadata_file = "RefData/LANL_Qry_Env_GE1000_Feb18_2025.xlsx"

aligned_seqs_output_file = "Aligned_Seqs_" + GENE + ".xlsx"
mutations_output_file_csv = "Mutations_" + GENE + ".csv"
mutations_output_file_excel = "Mutations_" + GENE + ".xlsx"

# Retrieve HXB2 gp120 or gp41 aas from GenBank. 
# gp120 has 511 aas; gp41 has 345 aas
hxb2_seq_path = "RefData/HXB2_" + GENE + "_AA.fasta"
(header, hxb2_seq) = read_fasta(hxb2_seq_path)
ref_seq_len = len(hxb2_seq)

# Retrive the AA seqs from the excel file with ~6000 filtered seqs
# File does not contain HXB2; # Headers are "Accession", "NumAAs", "Sequence"
df_seqs = pd.read_excel(env_sequences_file, engine='openpyxl')

# Retrieve the metadata from LANL with subtype from the 
# RefData/LANL_Qry_Env_GE1000_Feb18_2025.xlsx file
df_subtypes = pd.read_excel(metadata_file, sheet_name= 'Subtypes', engine='openpyxl')
df_seqs_w_subtypes = df_seqs.merge(df_subtypes[['Accession', 'Subtype']], on='Accession', how='left')
df_seqs_w_subtypes = df_seqs_w_subtypes[['Accession', 'NumAAs', 'Subtype', 'Sequence']]
dict_seqs = df_seqs_w_subtypes.set_index('Accession').to_dict(orient='index')

# Aligns sequences and writes them to a excel file that contains an index,
# the subtype, acccession number, numAAs, seq_len (query), score, 
# the aligned hxb2 sequence, andthe aligned query sequence
# This takes about 5 minutes 
def align_all_sequences_to_hxb2(dict_seqs):
    counter = 0
    alignment_output = []
    for acc, details in dict_seqs.items():
        subtype = details['Subtype']
        aas = details['Sequence']
        counter += 1
        if counter % 100 == 0:
            print(counter)
        (score, hxb2, query) = align_2_aa_seqs(hxb2_seq, aas)
        seq_len = len(query)
        if seq_len < (ref_seq_len - 50):
            continue
        alignment_output.append([counter, subtype, acc, seq_len, score, hxb2, query])
        columns = ["Counter", "Subtype", "Accession", "SeqLen", "Score", "HXB2", "Query"]
    alignment_df = pd.DataFrame(alignment_output, columns=columns)
    return alignment_df
    
alignment_df = align_all_sequences_to_hxb2(dict_seqs)
alignment_df.to_excel(aligned_seqs_output_file, index=False)
#sys.exit("Exiting the program.")

# Read from the aligned_sequences file to create a "mutation profile" spreadsheet 
# in which the amino acid at each position is compare to HXB2
df_aligned_seqs = pd.read_excel(aligned_seqs_output_file)
   
header_positions = [f"P{i}" for i in range(1,(ref_seq_len + 1))]
header_positions_string = ','.join(header_positions)
header1 = 'Acc,Subtype,Matches,' + header_positions_string + "\n"
header2 = 'K03455,B,' + str(ref_seq_len) + ',' + ','.join(hxb2_seq) + "\n"

counter = 0
with open(mutations_output_file_csv, 'a') as file:
    file.write(header1)
    file.write(header2)
    for index, row in df_aligned_seqs.iterrows():
        counter += 1
        acc = row['Accession']
        subtype = str(row['Subtype'])
        if subtype == 'nan':
            subtype = 'Unk'
        seq_len = row['SeqLen']
        aligned_hxb2 = row['HXB2']
        aligned_query = row['Query']
    
        num_missing_pos = check_hxb2_start(aligned_hxb2, ref_seq_len)
        if num_missing_pos !=0:
            leader_aas = '.' * num_missing_pos
            aligned_hxb2 = leader_aas + aligned_hxb2
            aligned_query = leader_aas + aligned_query
        query_mut_dict = compare_2_strings(hxb2_seq, aligned_hxb2, aligned_query)
        num_matches = str(sum(1 for value in query_mut_dict.values() if value == '-'))

        #print(query_mut_dict)
        values = [acc, subtype, num_matches]
        for key in range(1, (ref_seq_len + 1)):  
            if key in query_mut_dict:  # Check if the key exists in the dictionary
                values.append(query_mut_dict[key])
        result_string = ','.join(values) + "\n"
        #print(f"{counter},{result_string}")
        file.write(result_string)

    # Convert the csv file to an excel file
    df = pd.read_csv(mutations_output_file_csv)
    df.to_excel(mutations_output_file_excel, index=False, engine='openpyxl')

    # Delete the CSV file
    if os.path.exists(mutations_output_file_csv):
        os.remove(mutations_output_file_csv)
        print(f"CSV file '{mutations_output_file_csv}' has been deleted.")
    else:
        print(f"CSV file '{mutations_output_file_csv}' does not exist.")
