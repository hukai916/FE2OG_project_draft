"""
To select out the 50 target clusters and rename the clusters.
python rename_select.py /Users/Kai/Desktop/GitHub/FE2OG_project_draft/Blastp_cluster/filter_0.95_0.35/114_50.csv
"""

import sys
import subprocess
from pathlib import Path

filename = sys.argv[1]
work_dir = Path(__file__).resolve().parent.parent

input_dir = str(work_dir) + "/Blastp_cluster/filter_0.95_0.35/Cluster/"
output_dir = str(work_dir) + "/Blastp_cluster/filter_0.95_0.35/Select50/"

for line in open(filename):
    tem = line.split(",")
    in_name = tem[0] + ".fasta"
    out_name = "Cluster_" + tem[2] + '_' + tem[1] + '_' + tem[-1][:-1] + ".fasta"
    command = "cp " + input_dir + in_name + " " + output_dir + out_name
    subprocess.call(command, shell=True)

print(output_dir)
