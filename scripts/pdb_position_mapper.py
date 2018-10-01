"""
To map the TMalign coordinates onto the site positions from MSA.
Usage:
python pdb_position_mapper.py
"""
from subprocess import Popen, PIPE, STDOUT, run
from Bio import SeqIO


def site_mapper(pdb_seq, msa_seq, siteList):
# POPEN.run is an alternative to Popen() and P.communicate()
# This function will return 1-based coordiates in MSA based on the positions in pdb seq
    input_template = ">1\n" + pdb_seq + '\n' + ">2\n" + msa_seq + '\n'
    proc = run(['clustalo', '--infile=-'], stdout=PIPE, input=input_template, stderr=STDOUT, encoding='ascii')
    temList = proc.stdout.split("\n")
    #print(temList)
    template = ''.join(temList[1: int(len(temList)/2)])
    target   = ''.join(temList[int(len(temList)/2 + 1):])
    #template, target = proc.stdout.split('\n')[1], proc.stdout.split('\n')[3]
    alignment_indicator = ['WRONG!' if "-" in target else '']
    count = 0
    siteListMap = []
    if alignment_indicator[0]: siteListMap.append(alignment_indicator[0])
    for i in range(len(template)):
        if not template[i] == '-':
            count = count + 1
            if count in siteList: siteListMap.append(i + 1)
    print("site_mapper:\n", temList)
    return(siteListMap)

def chimera_mapper(chimerafilename, siteList): # Given siteList (fas-based VioC_18) and Chimera alignment, output fas-based target sites.
    seqList     = [record.seq for record in SeqIO.parse(chimerafilename, "fasta")]
    VioC_seq    = str(seqList[0]).replace('.', '-') # sometimes, Chimera will use dot instead of - for the gapped pair alignment.
    Target_seq  = str(seqList[1]).replace('.', '-')
    #print(Target_seq)
    siteLocator = 0
    temList     = []
    for i in range(0, len(VioC_seq)):
        if not VioC_seq[i] == '-':
            siteLocator = siteLocator + 1
            if siteLocator in siteList:
                temList.append(i+1)
    resList = []
    posLocator = 0
    for i in range(0, len(Target_seq)):
        if i in temList:
            resList.append(posLocator)
        if not Target_seq[i] == '-': posLocator = posLocator + 1
        #print(i, Target_seq[i])

    #print(resList, "???????")
    #print(temList)
    #print(siteList)
    #print(resList)
    return(resList, temList)

def chimera_mapper2(chimerafilename, currChi):
    seqList     = [record.seq for record in SeqIO.parse(chimerafilename, "fasta")]
    VioC_seq    = str(seqList[0]).replace(".", '-')
    Target_seq  = str(seqList[1]).replace(".", '-')
    resList = []
    posLocator = 0
    #print(currChi, "???")
    #print(type(currChi))
    #print(type(currChi.split()))
    currChi = [int(x) if not x == 'gap' else x for x in currChi.split()]
    #print(currChi)
    #currChi = [int(x) if not x == 'gap' for x in currChi.split()]
    for i in range(0, len(Target_seq)):
        if i in currChi:
            resList.append(posLocator)
        if not Target_seq[i] == '-': posLocator = posLocator + 1
    #print(resList)
    return(resList)


def msa_mapper(mapped_sites, msa_file): # This is to retrieve the position in the MSA based on the coordinates in mapped_sites
    msa_sites = []
    for record in SeqIO.parse(msa_file, 'fasta'):
        count = 0
        for i in range(0, len(record.seq)):
            if not record.seq[i] == '-':
                count = count + 1
                if count in mapped_sites: msa_sites.append(i+1)
        break
    return(msa_sites)


if __name__ == "__main__":
    site_mapper('TAGATTT', 'TTTAACGTTT', [1,2,4])
