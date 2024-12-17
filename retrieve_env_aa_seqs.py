import csv
import re
import sys

from hiv_seq_utils import (read_fasta, safe_fetch_env_aa)

# Create a CSV file with the gp160 amino acid sequences of all the entries in 
# the LANL filtered envelope dataset (about 6500 accession numbers).
# This has to be run just once. The program takes several hours because GenBank
# responds slowly.
# The file will not contain HXB2 because that will be used as the reference
# sequence in Env_align.py

# Of note, many AA sequences are unable to be returned for various reasons 
# Even HXB2 AAs were not returned (K03455). 
# Further optimization of the program is required to better control what AAs
# are identified. Translation of NAs may be required in some instances. In others,
# it may be useful to return why AAs were not present  
# Overall, 5862 sequences of lengths 462 to 901 were returned. These included
# 50 sequences between 462 and 809 AAs.
# 750 were missing; 44 were shorter than 450 AAs; These were deleted from the 
# the excel file but left in the CSV file.


lanl_filtered_acc = "RefData/LANL_HIV1_FLT_2021_accessions.txt"
acc_numbers = []
with open(lanl_filtered_acc, 'r') as file:
    for line in file:
        acc_number = line.strip()
        acc_numbers.append(acc_number)

# Get the envelope amino acid sequences
gp160_aa_sequences = []
counter=0
for acc in acc_numbers:
    counter += 1
    if counter % 100 == 0:
        print(counter)
    (protein_id, aas) = safe_fetch_env_aa(acc)
    if  protein_id is not None and aas is not None:
        gp160_aa_sequences.append(f"{acc},{len(aas)},{aas}\n")
    else:
        gp160_aa_sequences.append(f"{acc},-,-\n" )
        print(acc)

# Write the entries to a text file
full_text = ''.join(gp160_aa_sequences)
with open('gp160_aa_sequences.csv', 'w') as file:
    file.write(full_text)







