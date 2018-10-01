"""
This script is to parse the .clstr file together with cluster id file
in order to reorganize the order of sequences in each cluster by putting
the structure-associated seqs on top of each cluster.
Usage:
python parse_clstr.py filter20.clstr Fatcat/153_0.75_domain.txt
"""
from Bio import SeqIO
import sys
import re
import os

clstr_file = sys.argv[1]
cluster_file = sys.argv[2]
all_seq_file = sys.argv[3]
pdb_list = [x.split()[0].upper() for x in open(cluster_file)]
pdb_list = [x[:-1] + "_" + x[-1] for x in pdb_list]
all_dict = {}
for record in SeqIO.parse(all_seq_file, 'fasta'):
    key = record.description.split()[0]
    all_dict[key] = []
    all_dict[key].append(record.description)
    all_dict[key].append(record.seq)

def reorder(key, seq_list):
    seqList = [x[1:-3] for x in seq_list]
    check = sorted(set(seqList).intersection(set(pdb_list)))
    reorder_name = []
    if check: # some clusters don't contain sequences that are associated with pdb structures. 92 clusters left.
        reorder_name = check
        for x in seqList:
            if not x in reorder_name:
                reorder_name.append(x)
    if len(reorder_name) > 0:
        return(reorder_name)

clstrDict = {}

for line in open(clstr_file):
    if line.startswith(">"):
        key = '_'.join(line.split())
        clstrDict[key] = []
    else:
        clstrDict[key].append(line.split()[2])

for key in clstrDict:
    seq_list = clstrDict[key]
    res_reorder = reorder(key, seq_list)
    if res_reorder:
        tem = '/'.join(clstr_file.split("/")[:-1])
        _filename = tem + '/Cluster/' + key[1:] + '.fasta'
        if not os.path.exists(os.path.dirname(_filename)): os.makedirs(os.path.dirname(_filename))
        outFile = open(_filename, 'w+')
        print(tem + '/' + key[1:] + '.fasta')
        for seqName in res_reorder:
            outName = '>' + all_dict[seqName][0]
            outSeq  = all_dict[seqName][1]
            outFile.write(outName + '\n' + str(outSeq) + '\n')
        outFile.close()

        #break
    #for x in seq_list: print(x)
