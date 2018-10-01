"""
A Python wrapper for cd-hit: a fast sequence clustering package.
"""

import sys
import subprocess
from pathlib import Path

filename = sys.argv[1]
cutoff   = sys.argv[2] # input the clustering cutoff
work_dir = Path(__file__).resolve().parent.parent

queryName = filename.split("/")[-1].split(".")[0]
inFile   = str(work_dir) + '/Blastp_fasta/' + filename.split("/")[-1]
outFile  = str(work_dir) + "/Blastp_cluster/" + queryName + '.cluster'

cdhit_command = ' '.join(['cd-hit', '-i', inFile, '-o', outFile, '-c', str(cutoff), '-n 5 -M 1000 -d 0 -T 2'])

proc_cdhit = subprocess.Popen(cdhit_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
out, err = proc_cdhit.communicate()
out = out.decode().split("\n")
totalHit = out[-9].split()[0]
totalCluster = out[-9].split()[2]

logfile = open(''.join([str(work_dir), '/Blastp_cluster/', 'log.txt']), 'a+')
loginfo = ' '.join([queryName, totalHit, totalCluster])
print(loginfo) # queryName totalHit_num totalCluster_num
logfile.write(loginfo + "\n")
