"""
This wrapper will compile the scripts used in the entire pipeline.
Usage:
python wrapper.py ***.rename.clustal
The output will be saved in MSA and PDF folder.
Note, tm_pos is the coordinate based on tm-align result, basically is the structurally solved sequences.
Note, pdb_pos is the coordinate displayed in Chimera structure viewer.
Note, fa_pos is the coordinate based on fasta sequence (the first one in the MSA).
tm_pos = pdb_pos - start_pos(excluding non-solved seq) + 1
"""
from Bio import SeqIO
from tmalign import tm_align_pair_wrapper, chimera_align_pair_wrapper
from pdb_position_mapper import site_mapper, msa_mapper
import subprocess # Note that run has substituted call
import sys
import os

filename = sys.argv[1]

clusterID, pdbID = filename.split("_")[-3:-1]
for record in SeqIO.parse(filename, "fasta"): # Get the seq in MSA
    msa_seq = ''.join(str(record.seq).split("-"))
    break

res = tm_align_pair_wrapper(target=pdbID)[1]
print("from TMalign:", res)

res = chimera_align_pair_wrapper(target=pdbID)
print("from Chimera:", res)


# return value from tm_align_pair_wrapper: [query_seq_withoutgap, ref_seq_withoutgap, query_id, ref_id, sites_aligned, sites_target_pos]

pdb_seq         = res[0] # Get the seq from pdb file, only the structually solved part will be in the seq
ref_seq         = res[1] # Get the seq from ref (template) pdb file, only the structurally solved part will be in the seq
siteTargetList  = res[-1]# Get the site coordinates for target seq from the TM-aligned output(structurally solved part)
siteAlignInfo   = res[-2]# Get the structurally alignment info for the selected site pairs

# Note, if TM-align didn't output the best equivalent sites, we still need to manually adjust the site list before continuing the following steps.
# To manually adjust the siteTargetList, get clues from Chimera structure alignment, the pdb_pos can be converted to tm_pos by: pdb_pos - start_pos(excluding non-solved) + 1 = tm_pos.

#print(ref_seq)
#print(pdb_seq)

mapped_site = site_mapper(pdb_seq, msa_seq, siteTargetList) # get the mapped site coordinates (1-based) in msa_seq
#msa_site    = msa_mapper(mapped_site) # to retrieve 1-based coordinates in msa file
msa_site    = msa_mapper(mapped_site, filename)

msa_list = []

for record in SeqIO.parse(filename, 'fasta'):
    temMSA = ''.join(record.seq[i] if not record.seq[i] == '-' else 'X' for i in msa_site)
        # replace gap, '-', with 'X'
    msa_list.append(temMSA)
print("\n")
print("Retrieve equivalent sites based on TMaligned structures:")
print('TM-alignment for each site pair: ":" is perfect match, "." is somewhat match\n', siteAlignInfo)

print('Target sites (tm_pos): TMalign result based coordinates\n', siteTargetList) # equivalent sites in target from TM_pos
print('Target sites (fa_pos): fasta seq based coordinates\n', mapped_site) # equivalent sites in target from fasta coordinates
print('Target sites (msa_pos): msa based coordinates\n', msa_site, "\n") # from msa coordinates

homedir = os.path.dirname(os.getcwd())

## Save MSA into MSA folder
msadir  = os.path.join(homedir, 'MSA')
msafilename = os.path.split(filename)[-1].split(".")[0] + '.msa.txt'
msafile0 = os.path.join(msadir, msafilename)
os.makedirs(os.path.dirname(msafile0), exist_ok=True)
msafile = open(msafile0, 'w+')
msafile.write('\n'.join(msa for msa in msa_list))
msafile.close()
print("Done retrieving equivalent MSA: saved in MSA folder.")

## Build Seqlogo using MSA produced above
logodir = os.path.join(homedir, 'LOGO')
logofilename = os.path.split(filename)[-1].split(".")[0] + '.logo'
logofile0 = os.path.join(logodir, logofilename)
os.makedirs(os.path.dirname(logofile0), exist_ok=True)
width = len(msa_site)
command = " ".join(["seqlogo -F PDF", "-f", msafile0, "-l 1 -m", str(width), "-o", logofile0, "-c -S -n -w", str(width+2)])
subprocess.run(command, shell=True)
print("Done generating letter LOGO: saved in LOGO folder.")


#seqlogo -F PDF -f MSA_18/SyrB2.msa -l 1 -m 18 -o PDF/sdg -c -S -n -w 20
#print([x-8-1 for x in mapped_site]) # to get the coordinate numbers in pdb file. plus starting pos then minus one.
