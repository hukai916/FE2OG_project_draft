"""
To rename the select 50 cluster queries.
Usage:
python XXX.fasta
"""

import sys
from pathlib import Path
from Bio import SeqIO

work_dir = Path(__file__).resolve().parent.parent

infile = sys.argv[1]
outfile = open(infile + ".rename", 'w+')

count = 1
seq_name = infile.split("/")[-1][8:-6]

for record in SeqIO.parse(infile, 'fasta'):
    seqName = ">C" + seq_name + "_" + str(count)
    count = count + 1
    outfile.write(seqName + "\n")
    outfile.write(str(record.seq) + "\n")

outfile.close()

#print(outfile)
print(seq_name)
