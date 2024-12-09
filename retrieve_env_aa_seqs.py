import csv
import re

from hiv_seq_utils import (read_fasta, safe_fetch_env_aa)

# Create a CSV file with the gp160 amino acid sequences of all the entries in 
# the LANL filtered envelope dataset. 

# The first row in the file will have the gp120 or gp41 amino acid HXB2 sequence
gp120_aa_sequences = []
gp41_aa_sequences = []
hxb2_gp120 = "RefData/HXB2_GP120_AA.fasta"
hxb2_gp41 = "RefData/HXB2_GP41_AA.fasta"
hxb2_na_acc = "K03455"
(gp120_header, gp120_sequence) = read_fasta(hxb2_gp120)
(gp41_header, gp41_sequence) = read_fasta(hxb2_gp41)
hxb2_gp120_entry = f"{hxb2_na_acc}, {len(gp120_sequence)}, {gp120_sequence}"
hxb2_gp41_entry = f"{hxb2_na_acc}, {len(gp41_sequence)}, {gp41_sequence}"

# LANL filtered viruses - text file with GenBank NA accessions for about 7000
lanl_filtered_acc = "RefData/LANL_HIV1_FLT_2021_accessions.txt"
acc_numbers = []
with open(lanl_filtered_acc, 'r') as file:
    for line in file:
        acc_number = line.strip()
        acc_numbers.append(acc_number)

# Get the envelope amino acid sequences
gp120_aa_sequences.append(hxb2_gp120_entry)
counter=0
for acc in acc_numbers:
    counter += 1
    if counter % 100 == 0:
        print(counter)
    (protein_id, aas) = safe_fetch_env_aa(acc)
    if  protein_id is not None and aas is not None:
        gp120_aa_sequences.append(f"{acc},{len(aas)},{aas}\n")
    else:
        gp120_aa_sequences.append(f"{acc},-,-\n" )
        print(acc)

# Write the entries to a text file
full_text = ''.join(gp120_aa_sequences)
with open('gp120_aa_sequences.csv', 'w') as file:
    file.write(full_text)







