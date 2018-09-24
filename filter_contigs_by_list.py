#!/usr/bin/python3
import sys
from Bio import SeqIO

fasta = sys.argv[1]
fasta_list = sys.argv[2]
corrected_file = sys.argv[3]

contigs_to_keep = []

for line in open(fasta_list):
    contig=line.rstrip()
    contigs_to_keep.append(contig)

with open(corrected_file, 'w') as corrected:
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if seq_record.id in contigs_to_keep:
            SeqIO.write(seq_record, corrected, 'fasta')
    
