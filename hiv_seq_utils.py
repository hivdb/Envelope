import csv
import re
from io import StringIO
from Bio import Align, Entrez, SeqIO
from Bio.Align import substitution_matrices

Entrez.email = "rshafer@stanford.edu"
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYBXZ")


def fetch_sequence(accession):
    #print("Accession:", accession)
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    gb_content = handle.read()
    handle.close()
    handle_gb = StringIO(gb_content)
    record = SeqIO.read(handle_gb, "genbank")
    handle_gb.close()
    fasta_seq = record.seq
    return gb_content, fasta_seq


def fetch_env_aa(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get('gene', ['Unknown gene'])[0]
            product = feature.qualifiers.get('product', ['Unknown gene'])[0]
            if gene_name == "env" or product[:3].lower == "env":
                protein_id = feature.qualifiers.get('protein_id', ['No protein ID'])[0]
                if protein_id != 'No protein ID':
                    protein_sequence = feature.qualifiers.get('translation', ['No protein translation'])[0]
                    #print(f"Gene: {gene_name}, ProteinID: {protein_id} Sequence: {protein_sequence}")
                    return protein_id, protein_sequence       


def safe_fetch_env_aa(acc):
    result = fetch_env_aa(acc)
    if result is None:
        return (None, [])
    return result


def read_fasta(file_path):
    # Open the FASTA file
    with open(file_path, 'r') as file:
        # Initialize an empty string to store the sequence
        sequence = ''
        # Iterate over each line in the file
        for line in file:
            # Strip newline characters from the line
            line = line.strip()
            # Check if the line is a header (which starts with '>')
            if line.startswith('>'):
                # Optionally, you can print or store the header if needed
                header = line[1:]  # This removes the '>' character
                #print(f"Header: {header}")
            else:
                # Append the line to the sequence
                sequence += line
        return header, sequence


def align_2_aa_seqs(seq1, seq2):
    seq2 = ''.join(seq2.split())
    invalid_chars1 = set(seq1) - valid_amino_acids
    invalid_chars2 = set(seq2) - valid_amino_acids
    if invalid_chars1:
        print("Invalid characters in HXB2:", invalid_chars1)
    if invalid_chars2:
        print("Invalid characters in Query:", invalid_chars2)
    aligner = Align.PairwiseAligner()
    matrix = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = matrix
    aligner.mode = 'local'  # Choose 'global' or 'local'
    aligner.open_gap_score = -8  # Penalty for opening a gap
    aligner.extend_gap_score = -1  # Penalty for extending a gap  
    score = aligner.score(seq1, seq2)
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]

    alignment_string = format(best_alignment)
    lines = alignment_string.split('\n')
     
    aligned_target = ""
    aligned_query = ""
    for i, line in enumerate(lines):
        if i % 4 == 0:  # Target lines occur every 4 lines starting from 0
            line = remove_numbers_from_end(line)  
            aligned_target += line.split()[-1]  # Assume last part after space is sequence
        if (i - 2) % 4 == 0:
            line = remove_numbers_from_end(line)
            aligned_query += line.split()[-1]  # Assume last part after space is sequence
    
    # Approximately 80 of 5000 sequences had trailing "target" or "query" strings
    aligned_target = remove_trailing_string(aligned_target, "target")
    aligned_target = remove_numbers_from_end(aligned_target)
    aligned_query = remove_trailing_string(aligned_query, "query")
    aligned_query = remove_numbers_from_end(aligned_query)
    return (best_alignment, score, aligned_target, aligned_query) 

def remove_trailing_string(sequence, string):
    if sequence.endswith(string):
        return sequence[:-len(string)]
    return sequence

def remove_numbers_from_end(s):
    # Pattern to match numbers at the end of the string
    return re.sub(r'\d+$', '', s)  


def compare_2_strings (hxb2_seq, ref, query):
    query_mut_dict = {}
    pos = 1
    for i, j in zip(ref, query):
        #print("Pos:", pos, " i:", i, " j:", j)
        if i == j:
            query_mut_dict[pos] = '-'
            pos += 1
        elif i != '-' and i != j and j != '-':
            query_mut_dict[pos] = j
            pos += 1
        elif j == '-':
            query_mut_dict[pos] = '.'
            pos += 1
        elif i == '-':
            if (query_mut_dict[pos-1] == "-"):
                query_mut_dict[pos-1] = hxb2_seq[pos-1]
            query_mut_dict[pos-1] = query_mut_dict[pos-1] + j
    return(query_mut_dict)      


def check_hxb2_start(aligned_hxb2, hxb2_num_aas):
    aligned_hxb2_wout_insertions = aligned_hxb2.replace('-','')
    num_missing_pos = hxb2_num_aas - len(aligned_hxb2_wout_insertions)
    return num_missing_pos