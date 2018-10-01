"""
This is the refined semi-automatic pipeline which
takes into account of the Chimera adjustment as well.
Usage:
python wrapper_refine.py XXX.clustal
Note that, chimera seq is the same as fas-based seq, which represent all sites in pdb file including unsolved sites.
"""

"""
1. Take two PDBs and TMalign them.
2.
"""

from Bio import SeqIO
from tmalign import tm_align_pair_wrapper, tm_align_pair
from Constant import VioC_18 as VioC_18_fas # selected sites in VioC based on fas coordinates.
from pdb_position_mapper import chimera_mapper, chimera_mapper2, msa_mapper
import sys
import os
import subprocess

filename = sys.argv[1]
homedir  = os.path.dirname(os.getcwd())
clusterID, pdbID = filename.split("_")[-3:-1]


tmdir    = os.path.join(homedir, 'TMres')
tmfilename = os.path.split(filename)[-1].split(".")[0] + '.tm.fasta'
tmfile0 = os.path.join(tmdir, tmfilename)
os.makedirs(os.path.dirname(tmfile0), exist_ok=True)

resAll = tm_align_pair_wrapper(target=pdbID)

res = resAll[1]
#print(res[3][0]) # TMalign -- query
#print(res[3][1]) # TMalign -- target

tmfile = open(tmfile0, 'w+')
tmfile.write('\n'.join(['>VioC', res[3][0], '>' + pdbID, res[3][1]]))
tmfile.close()

#print("from TMalign:", res)

# Below to examine chimera output alignment.
chimerafilename = os.path.split(filename)[-1].split(".")[0] + '.tm.chimera.fasta'
chimerafilename = os.path.join(tmdir, chimerafilename)

def MSA_LOGO(msa_pos, currPara = "default"): # retrieve sites for msa file, save in MSA folder and LOGO in LOGO folder. curr is to indicate the positions of 'gap'
    ## Save MSA into MSA folder
    msadir  = os.path.join(homedir, 'MSA')
    msafilename = os.path.split(filename)[-1].split(".")[0] + '.msa.txt'
    msafile0 = os.path.join(msadir, msafilename)
    os.makedirs(os.path.dirname(msafile0), exist_ok=True)

    msa_list = []
    print(msa_pos)
    for record in SeqIO.parse(filename, 'fasta'):
        #print(record.seq)
        #print(msa_pos)
        temMSA = ''.join(record.seq[i-1] if not record.seq[i-1] == '-' else 'X' for i in msa_pos)

        # replace gap, '-', with 'X'
        msa_list.append(temMSA)

# Modify the msa by adding in the X for gapped positions.
    msafile = open(msafile0, 'w+')
    gapPos = [i for i in range(0, len(currPara.split())) if currPara.split()[i] == 'gap']
    for i in range(0, len(msa_list)):
        for gap in gapPos:
            msa_list[i] = msa_list[i][0:gap] + "X" + msa_list[i][gap:]
    print("This is gap position:    ", currPara)
    print("This is msa_pos:    ", msa_pos)
    msafile.write('\n'.join(msa for msa in msa_list))
    msafile.close()
    print("Done retrieving equivalent MSA: saved in MSA folder.")

    ## Build Seqlogo using MSA produced above
    logodir = os.path.join(homedir, 'LOGO')
    logofilename = os.path.split(filename)[-1].split(".")[0] + '.logo'
    logofile0 = os.path.join(logodir, logofilename)
    os.makedirs(os.path.dirname(logofile0), exist_ok=True)
    width = max(len(msa_pos), len(currPara.split()))
    command = " ".join(["seqlogo -F PDF", "-f", msafile0, "-l 1 -m", str(width), "-o", logofile0, "-c -S -n -w", str(width+2)])
    #print(command)
    subprocess.run(command, shell=True)
    print("Done generating letter LOGO: saved in LOGO folder.")

def chimeraPos(chimera_pos): # Save chimera_pos into files.
    chidir  = os.path.join(homedir, 'Chimera_pos')
    chifilename = os.path.split(filename)[-1].split(".")[0] + '.chimera.txt'
    chifile0 = os.path.join(chidir, chifilename)
    os.makedirs(os.path.dirname(chifile0), exist_ok=True)

    chifile = open(chifile0, 'a+')
    chifile.write(' '.join(str(chi) for chi in chimera_pos))
    chifile.write("\n")
    chifile.close()

if os.path.exists(chimerafilename):
    #print("VioC_18_mas:\n", VioC_18_fas)
    #print("from tmaling:\n", resAll[0])
    fas_pos, chimera_pos  = chimera_mapper(chimerafilename, VioC_18_fas)
    #print("?????")
    #print(chimerafilename)
    #print(fas_pos)
    #print(chimera_pos)
    #print(chimerafilename)
    msa_pos  = msa_mapper(fas_pos, filename)
    print("VioC_site (fas_pos):\n", *VioC_18_fas)
    print("Chimera_pos:\n", *chimera_pos)
    print(pdbID + "_site (fas_pos):\n", *fas_pos)
    print(pdbID + "_site (msa_pos):\n", *msa_pos)
    chimeraPos(chimera_pos) # create a file that contains chimera positions.
    print("\nBelow is TM-alignment aided MSA_LOGO():\n")
    #tem_msa = " ".join(map(str, msa_pos))
    MSA_LOGO(msa_pos)

else:
    print("NO! Not provided with Chimera file yet! Provide below:\n", chimerafilename, "\n")

# Below is to use adjusted Chimera pos file. Note, the updated coordinates should be at the second line of the chifile.
chidir  = os.path.join(homedir, 'Chimera_pos')
chifilename = os.path.split(filename)[-1].split(".")[0] + '.chimera.txt'
chifile0 = os.path.join(chidir, chifilename)

if os.path.exists(chifile0):
    chiposList = [line.split("\n")[0] for line in open(chifile0)]
    curr       = [chiposList[0] if len(chiposList) <2 else chiposList[1]][0]
    original   = chiposList[0]
    print("Original Chi_pos:\n", original)
    print("Current Chi_pos:\n", curr)
    if not curr == original:
        gapPos = [i for i in range(0, len(curr.split())) if curr.split()[i] == 'gap']
        fas_pos = chimera_mapper2(chimerafilename, curr)
        fas_pos_print = fas_pos
        for gap in gapPos: fas_pos_print = fas_pos_print[0:gap] + ["gap"] + fas_pos_print[gap:]

        print("fas_pos:\n", *fas_pos_print)
        msa_pos = msa_mapper(fas_pos, filename)
        msa_pos_print = msa_pos
        for gap in gapPos: msa_pos_print = msa_pos_print[0:gap] + ["gap"] + msa_pos_print[gap:]
        print("mas_pos:\n", *msa_pos_print)
        print("\nBelow is Chimera-aided MSA_LOGO(): \n")

        MSA_LOGO(msa_pos, currPara=curr)

        print("Done with updated Chimera pos!")
    else:
        print("OK! No updated chimera_pos file detected! Provide the second line in the following file if any change to the first line:\n Chimera_pos/Cluster_X_XXXXXX.chimera.txt")

else:
    print("OK! No updated chimera_pos file detected! Provide the second line if necessary:\n Chimera_pos/Cluster_X_XXXXXX.chimera.txt")
