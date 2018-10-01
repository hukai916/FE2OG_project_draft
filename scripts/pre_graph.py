"""
Process raw_pdb_comparison.txt file by using smaller TMscore for each pair.
Usage:
python pre_graph.py ../PDB/raw_pdb_comparison.txt > ../PDB/tm_score.txt
"""

import sys
filename = sys.argv[1]

minDict = {}

# Find the minimum tmScore value for each pair.
for line in open(filename):
    if len(line.split()) == 3:
        key1 = '+'.join([line.split()[0], line.split()[1]])
        key2 = '+'.join([line.split()[1], line.split()[0]])
        if not (key1 in minDict or key2 in minDict):
            minDict[key1] = []
            minDict[key1].append(float(line.split()[2]))
        else:
            if key1 in minDict:
                minDict[key1].append(float(line.split()[2]))
            else: minDict[key2].append(float(line.split()[2]))

# Assign the minimum value to each original pair.
for line in open(filename):
    if len(line.split()) == 3:
        key1 = '+'.join([line.split()[0], line.split()[1]])
        key2 = '+'.join([line.split()[1], line.split()[0]])
        if key1 in minDict:
            value = minDict[key1]
        if key2 in minDict:
            value = minDict[key2]
        print(*line.split()[0:2], min(value))


#for key in minDict:
#    print(key, minDict[key])

"""
for key in minDict:
    if len(minDict[key]) > 1: # get rid of self comparison
        if min(minDict[key]) > 0.5:
            print(*key.split("+"), min(minDict[key]))
    #print((minDict[key]))
"""
