"""
Read in Blastp output .xml file and retrieve a list of
fasta file sequences that satisfy:
sequence length between 0.5*len(matched_query) - 1.5*len(matched_query) and identity to query > 35%.
Usage:
python parseBlastp.py xxx.xml
output will be saved in Blastp_fasta folder.
This is a faster version.
"""

from Bio.Blast import NCBIXML
import sys
import subprocess
from pathlib import Path

def parseID(filename):
    result = open(filename, "r")
    records= NCBIXML.parse(result)
    idList = []
    for item in records:
        query_length = float(item.query_letters)
        for alignment in item.alignments:
            hit_length = float(alignment.length)
            #if 0.5 < hit_length / query_length < 2.0: # this filter is not valid since a lot of hits are imbedded in much larger proteins and would be filtered out according to this cutoff. Need to add in the hit_start and end.
            #hit_name   = alignment.id + alignment.def
            #print(dir(alignment))
            hit_def = alignment.hit_def.split(">")[0]
            for hsp in alignment.hsps:
                hit_start = hsp.sbjct_start
                hit_end   = hsp.sbjct_end
                if 0.5 < float(hit_end - hit_start) / query_length < 2.0:
                    if float(hsp.identities)/query_length > 0.35:
                        hit_id = alignment.accession
                        hit_seq = ''.join(letter if letter.isalpha() else '' for letter in hsp.sbjct)
                        hit_name = '>' + hit_id + ' ' + hit_def
                        #print(hit_name)
                        idList.append((hit_id, hit_start, hit_end, hit_length, hit_name, hit_seq))
    return(idList)

def formatFasta(idList):
    fastaList = []
    for id_tuple in idList:
        id = id_tuple[0]
        hit_start = id_tuple[1]
        hit_end   = id_tuple[2]
        hit_length= id_tuple[3]
        hit_name  = id_tuple[4]
        hit_seq   = id_tuple[5]
        nameline = hit_name + '\t' + str(int(hit_length)) + ':' + str(hit_start) + ':' + str(hit_end)
        seqline  = hit_seq
        fastaList.append(nameline + '\n' + seqline + '\n')
    return(fastaList)

def main():
    logfile = open(''.join([str(work_dir), '/Blastp_fasta/', 'log.txt']), 'a+')
    idList = parseID(filename)
    fastaList = formatFasta(idList)
    queryFile = str(work_dir) + '/Fasta/' + filename.split("/")[-1].split(".")[0] + ".fasta"
    queryList = '\n'.join([line.split("\n")[0] for line in open(queryFile)]) + '\n'
    fastaList = [queryList] + fastaList
    outfile = open(''.join([str(work_dir), '/Blastp_fasta/', filename.split("/")[-1], '.fasta']), 'w+')
    outfile.write(''.join(fastaList))
    outfile.close()
    loginfo = (" ".join([filename, ":", str(len(fastaList)), 'filtered blastp hits saved to Blastp_fasta folder. Query has been put on top.']))
    print(loginfo)
    logfile.write(loginfo + '\n')
    
if __name__ == "__main__":
    filename = sys.argv[1]
    work_dir = Path(__file__).resolve().parent.parent
    main()
