"""
Filter all_seq.fasta file, keep the unique ones based on the sequence and protein id.
Only one of the several varying-length proteins (the one with smallest length) with the same protein ID will be kept.
"""

import sys
from Bio import SeqIO

filename = sys.argv[1]
resDict  = {}


for record in SeqIO.parse(filename, 'fasta'):
    proteinID = record.description.split()[0]
    if not proteinID in resDict:
        resDict[proteinID] = [record.description, record.seq]
    else:
        if len(record.seq) < len(resDict[proteinID][1]):
            resDict[proteinID] = [record.description, record.seq]
for key in resDict:
    print('>' + resDict[key][0])
    print(resDict[key][1])
