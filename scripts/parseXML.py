import xml.etree.ElementTree as ET
import subprocess
import sys
from pathlib import Path
from multiprocessing import Pool

work_dir = Path(__file__).resolve().parent.parent

def xml2fasta(pdb): # retrieve xml file based on pdb, then parse out the fasta for the target chain.
    pdbid   = pdb[:-1].upper()
    chainid = pdb[-1].upper()
    xml_command = ''.join(['http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=', pdbid, '&customReportColumns=structureId,chainId,sequence&service=wsfile$format=xml'])

    proc_xml = subprocess.Popen(["curl", "-s", xml_command], stdout=subprocess.PIPE)
    out, err = proc_xml.communicate()
    tree = ET.ElementTree(ET.fromstring(out))
    root = tree.getroot()

    for record in root: # each record is a chain object
        chainid_r  = record.find('dimEntity.chainId').text.upper() # get the id of the chain
        pdbid_r    = record.find('dimEntity.structureId').text.upper()
        chainseq_r = record.find('dimEntity.sequence').text.upper()
        if chainid_r == chainid and pdbid_r == pdbid:
            #Optional, save retrieved fasta file into folder: Fasta
            filename = open(''.join([str(work_dir), '/Fasta/', pdbid_r, '_', chainid_r, '.fasta']), 'w+')
            filename.write('\n'.join(['>' + pdbid_r + '_' + chainid_r, chainseq_r]))
            filename.close()

            return(pdbid_r, chainid_r, chainseq_r)

def blastp(res_xml): # Blastp using fasta from res_xml
    pdbid    = res_xml[0]
    chainid  = res_xml[1]
    chainseq = res_xml[2]
    blastp_command = 'blastp -db nr -max_target_seqs 5000 -outfmt 5 -remote'
    proc_blastp = subprocess.Popen(blastp_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc_blastp.communicate(chainseq.encode()) # In Python2, simply proc_blastp.communicate(res[2]).
    #Optional, save retrieved blastp results into folder: Blastp
    filename = open(''.join([str(work_dir), '/Blastp/', pdbid, '_', chainid, '.xml']), 'w+')
    filename.write(out.decode())
    filename.close()

    return(out, err)

def main(pdb):
    #pdb = sys.argv[1]
    res_xml = xml2fasta(pdb)
    print(''.join([res_xml[0] + '_' + res_xml[1], '.fasta file saved to Fasta folder...']))
    print("Running BLASTP for " + res_xml[0] + "...")
    res_blastp = blastp(res_xml)
    print("BLASTP done for " + res_xml[0] + "!")

if __name__ == "__main__":
    filename = sys.argv[1]
    pdbList = [x.split()[0] for x in open(filename)]
    #print(pdbList)
    pool = Pool(1) # for parallel jobs, the parallel wont work for no reason!
    for x in zip(pdbList): print(x[0])
    pool.map(main,pdbList)
