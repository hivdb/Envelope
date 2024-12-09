import csv
import re
import sys
import pandas as pd

from hiv_seq_utils import (read_fasta, 
                           align_2_aa_seqs, 
                           compare_2_strings, 
                           check_hxb2_start)

pd.set_option('display.max_rows', 50)

# Retrieve HXB2 gp120 aas and its accession numbers from GenBank. It has 511 aas
hxb2_seq_path = "RefData/HXB2_GP120_AA.fasta"
hxb2_na_acc = "K03455"
hxb2_aa_acc = "AAB50262"
(header, hxb2_seq) = read_fasta(hxb2_seq_path)
hxb2_num_aas = len(hxb2_seq)
hxb2_fasta = f">{hxb2_na_acc}|{hxb2_aa_acc}\n{hxb2_seq}"

# Retrive the AA seqs from the excel file with 7000 filtered seqs
sequences_file = "RefData/LANL_FilteredAASeqs_june27.xlsx"
df = pd.read_excel(sequences_file, header=None, engine='openpyxl')
df.columns = ['NA_Acc', 'AA_Acc', 'NumAAs', 'Seq']
df = df[df['NumAAs'] > 600]
df_seqs = df[df['NA_Acc'] != 'K03455']

# Retrieve the metadata from LANL with subtype from the gp120_search_1pp (June 2022)
metadata_file = "RefData/LANL_gp120_search_1pp_13075_w_subtypes_June22.xlsx"
df_subtypes = pd.read_excel(metadata_file, sheet_name= 'Subtypes', engine='openpyxl')

df_seqs_w_subtypes = df_seqs.merge(df_subtypes[['NA_Acc', 'Subtype']], on='NA_Acc', how='left')
df_seqs_w_subtypes = df_seqs_w_subtypes.drop(['AA_Acc', 'NumAAs'], axis=1)
df_seqs_w_subtypes = df_seqs_w_subtypes[['NA_Acc', 'Subtype', 'Seq']]
dict_seqs = df_seqs_w_subtypes.set_index('NA_Acc').to_dict(orient='index')
#print(dict_seqs)

# Aligns sequences and writes them to a excel file that contains an index,
# the subtype, acccession number, seq_len, score, the aligned hxb2 sequence,
# the aligned query sequence
# This takes about 5 minutes and can be skipped to save time by commenting out
# the call to this function
def align_all_sequences_to_hxb2(dict_seqs):
    counter = 0
    alignment_output = []
    for acc, details in dict_seqs.items():
        subtype = details['Subtype']
        aas = details['Seq']
        counter += 1
        (alignment, score, hxb2, query) = align_2_aa_seqs(hxb2_seq, aas)
        seq_len = len(query)
        alignment_output.append([counter, subtype, acc, seq_len, score, hxb2, query])
        columns = ["Counter", "Subtype", "Accession", "SeqLen", "Score", "HXB2", "Query"]
    alignment_df = pd.DataFrame(alignment_output, columns=columns)
    return alignment_df
    
#alignment_df = align_all_sequences_to_hxb2(dict_seqs)
#alignment_df.to_excel("aligned_sequences.xlsx", index=False)
#sys.exit("Exiting the program.")


# The remaining code reads the data from the aligned_sequences file to create a
# "mutation profile" spreadsheet in which the amino acid at each position is compared
# to that of HXB2
aligned_seqs_file = "aligned_sequences.xlsx"
df_aligned_seqs = pd.read_excel(aligned_seqs_file)
   
header_positions = [f"P{i}" for i in range(1,512)]
header_positions_string = ','.join(header_positions)
header = 'Acc,Subtype,Matches,' + header_positions_string + "\n"
header2 = 'K03455,B,511,' + ','.join(hxb2_seq) + "\n"

counter = 0
with open("mutations.csv", 'a') as file:
    file.write(header)
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
    
        num_missing_pos = check_hxb2_start(aligned_hxb2, hxb2_num_aas)
        if num_missing_pos !=0:
            leader_aas = '.' * num_missing_pos
            aligned_hxb2 = leader_aas + aligned_hxb2
            aligned_query = leader_aas + aligned_query
        query_mut_dict = compare_2_strings(hxb2_seq, aligned_hxb2, aligned_query)
        num_matches = str(sum(1 for value in query_mut_dict.values() if value == '-'))

        print(query_mut_dict)
        values = [acc, subtype, num_matches]
        for key in range(1, 512):  
            if key in query_mut_dict:  # Check if the key exists in the dictionary
                values.append(query_mut_dict[key])
        result_string = ','.join(values) + "\n"
        print(f"{counter},{result_string}")
        file.write(result_string)


