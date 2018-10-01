"""
Align all to all structure alignment
"""
import subprocess
from subprocess import PIPE, STDOUT, Popen
import sys
import os
from pathlib import Path
from Constant import VioC_18, pdbList
from pdb_position_mapper import site_mapper
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

work_dir = Path(__file__).resolve().parent.parent
work_path= str(work_dir) + "/PDB/"

def tm_align_all(pdbList):
    #TMalign multiple structurs
    pairList = []
    homedir = os.path.dirname(os.getcwd())
    pdbFilePath = os.path.join(homedir, "PDB")

    for index, ref_pdb in enumerate(pdbList):
        for query_pdb in pdbList[:]:
            try:
                ref_pdb_t = os.path.join(pdbFilePath,ref_pdb)
                query_pdb_t = os.path.join(pdbFilePath, query_pdb)

                tm_output = subprocess.check_output(['TMalign', ref_pdb_t, query_pdb_t])
                #print(tm_output)
                res = parse_tm_output(tm_output, ref_pdb, query_pdb)
            except:
                continue
            pairList.append(res)
    return(pairList)

def tm_align_pair(pairList):
    #TMalign single pair of pdb structures, the first in the list will be reference
    try:
        tm_output = subprocess.check_output(['TMalign', pairList[0], pairList[1]])
        #print(tm_output)
        res = parse_tm_output(tm_output, pairList[0], pairList[1])
        #print("parse_tm_output:", res)
        return(res)
    except:
        pass

def chimera_align_pair(pairList):
    chimera_output = "/Users/Kai/Desktop/GitHub/FE2OG_project_draft/PDB/VioC_5LSQ_tmalign_match.fasta" # should be modifided to integrate Chimera automatically.
    res = parse_chimera_output(chimera_output, pairList[0], pairList[1])
    return(res)

def parse_chimera_output(chimera_output, ref_pdb, query_pdb):
    alignment = []
    for record in SeqIO.parse(chimera_output, 'fasta'):
        alignment.append(str(record.seq))
    alignment.append(':'*len(alignment[0]))
    return([ref_pdb, query_pdb, 1, alignment])





def parse_tm_output(tm_output, ref_pdb, query_pdb):
    #Parse tm_output from tm_align_all and extract tmscore and alignment
    #TMscore is normalizled to ref_pdb
    lines = tm_output.splitlines()
    lines = [x.decode("utf-8") for x in lines] # decode Bytes object to String

    #Extract tmscore
    tmscore = float()
    for line in lines:
        #print(type(line), line)
        if "TM-score=" in line:
            tmscore = float(line.split()[1])
            break
    #Extract alignment: in the order of ref-seq, query-seq, alignment
    alignment = [lines[-4], lines[-2], lines[-3]] #In the order of ref_seq, query_seq, alignment
    return([ref_pdb, query_pdb, tmscore, alignment])

def site_align(pair, ref_siteList):
    # pair is single result from tm_align_all output or output from tm_align_pair
    # ref_siteList refers to the actual site position (TM-aligned part) in in the fasta file, not the PDB coordinates
    alignment = pair[3][2]
    ref_seq   = pair[3][0]
    ref_seq_withoutgap = ''.join(pair[3][0].split("-"))
    query_seq = pair[3][1]
    query_seq_withoutgap = ''.join(pair[3][1].split("-"))
    ref_id    = pair[0]
    query_id  = pair[1]
    index = 0
    sites_aligned     = [] # Document the site pair aligment, : for execellent match while . for okay match
    sites_align_pos   = [] # Document the position of the sites using the coordinate from TM-aligned pair-wise alignment. 0 based nubmers.
    sites_target_pos  = [] # Document the position of the sites using the coordiante from target query. 0 based numbers.

    for i in range(len(ref_seq)):
        if not ref_seq[i] == '-':
            index = index + 1
            #print(ref_siteList, index, alignment[i])
            if index in ref_siteList: # sites in ref_siteList is 1-based.
                sites_aligned.append(alignment[i])
                sites_align_pos.append(i)
    index = 0
    for i in range(len(ref_seq)):
        if not query_seq[i] == '-':
            index = index + 1
        if i in sites_align_pos:
            sites_target_pos.append(index)

    #print(ref_seq)
    #print(query_seq)
    return([query_seq_withoutgap, ref_seq_withoutgap, query_id, ref_id, sites_align_pos, sites_aligned, sites_target_pos])
    #sum(1 for x in sites_align if x in (":", "."))

#pairList = ['VioC_VO.pdb', 'WelO5_5iqt.pdb']
#res = tm_align_pair(pairList)
#print(res)



def tm_align_pair_wrapper(template = 'VioC_VO.pdb', target = 'PDBID'): # the default value for template is VioC_VO.pdb
    pairList = [work_path + template, work_path + target + ".pdb"]
    res = tm_align_pair(pairList) # return tm-align result to res
    res2 = site_align(res, VioC_18) # further map tm-align result to find out MSA positions
    return(res2,res) # res2 contains: [query_seq_withoutgap, ref_seq_withoutgap, ref_id, query_id, sites_aligned, sites_target_pos]

def chimera_align_pair_wrapper(template = 'VioC_VO.pdb', target = 'PDBID'):
    pairList = [work_path + template, work_path + target + ".pdb"]
    res = chimera_align_pair(pairList) # return tm-align result to res
    res2 = site_align(res, VioC_18) # further map tm-align result to find out MSA positions
    return(res2) # res2 contains: [query_seq_withoutgap, ref_seq_withoutgap, ref_id, query_id, sites_aligned, sites_target_pos]

def site_main():
    pdbListFile = sys.argv[1]
    pdbListInput = [line.split()[0] for line in open(pdbListFile)]

    pdbSeq_path = str(work_dir) + "/Blastp_cluster/filter_0.95_0.35/Select50/Cluster*"

    for i in range(1,len(pdbListInput)):
        pairList = [work_path + 'VioC_VO.pdb', work_path + pdbListInput[i] + ".pdb"]
        res = tm_align_pair(pairList)
        target_name = pdbListInput[i]
        target_seq  = ''.join(res[3][1].split("-"))
        #site_mapper(target_seq, msa_seq, mapSite)

        #print(len(''.join(res[3][0].split("-"))))
        #pdprint(res[3][0])
        #print(pdbListInput[i], (''.join(res[3][1].split("-"))))
        #print(res[3][1])
        res2 = site_align(res, VioC_18)
        print(res2[1])
        """
        for j in range(1,len(pdbList)):
            pairList = [pdbList[i], pdbList[j]]
            res = tm_align_pair(pairList)
            res2 = site_align(res, VioC_18)
            print(res2)
        """
#site_main()

"""
res1 = tm_align_all(pdbList)
#print(res1[1])
print(res1[1][3][0])
print(res1[1][3][2])
print(res1[1][3][1])

        #print("\n".join([i for i in y[3]]))

res2 = [site_align(item, VioC_18, item[0]) for item in res1]
print(res2)
#print(res1)

"""
