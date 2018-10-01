"""
To download PDB files given the id list.
Usage:
python downloadPDB.py pdbid.list
"""

import subprocess
import sys
from pathlib import Path

work_dir = Path(__file__).resolve().parent.parent

filename = sys.argv[1]

for line in open(filename):
    pdb = line[:4]
    command = "wget " + "https://files.rcsb.org/download/" + pdb.upper() + ".pdb -O " + str(work_dir) + "/PDB/" + pdb.upper() + ".pdb"
    print(command)
    subprocess.call(command, shell=True)
