"""
Reorder the clusters by putting query seq on top.
"""
import sys
from pathlib import Path

fastaFile = sys.argv[1]
work_dir = Path(__file__).resolve().parent.parent

class ClstrOutput(object):
    def __init__(self, filename):
        self.parseClstr(filename)
    class cluster(object):
        def __init__(self, line):
            try:
                lineList = line.split()
                self.num = lineList[0]
                self.length = lineList[1][:-1]
                self.id  = lineList[2][1:-3]
                self.identity = lineList[3] if len(lineList) == 4 else lineList[4]
            except ValueError:
                msg = 'can not convert tokens {} to hit'
                raise ValueError(msg.format(line.split()))
    def parseClstr(self, filename):
        self.clusters = {}
        clusterCount  = -1
        for line in open(filename):
            tokens = line.split()
            if tokens[0] == ">Cluster": clusterCount = clusterCount + 1
            else:
                try: self.clusters[clusterCount].append(self.cluster(line))
                except KeyError: self.clusters[clusterCount] = [self.cluster(line)]

def getClusterRep(queryID, clstrFile):
    output = ClstrOutput(clstrFile)
    for clusterID in output.clusters:
        for seq in output.clusters[clusterID]:
            if seq.id == queryID:
                for seq in output.clusters[clusterID]:
                    if seq.identity == "*":
                        return(seq.id)
                        break
                break

def reorderFasta(repID, queryID, querySeq, clusterFile):
    outfile = ''.join([str(work_dir),  "/Blastp_cluster/", queryID, ".reorder.cluster"])
    temList = [line.split("\n")[0] for line in open(clusterFile)]
    resList = []
    for i in range(0, len(temList), 2):
        if not temList[i] == '>' + repID:
            resList.append(temList[i])
            resList.append(temList[i+1])
    resList = [''.join(['>', queryID]), *querySeq] + resList
    out = open(outfile, 'w+')
    out.write('\n'.join(resList))
    out.write('\n')
    out.close()

def main():
    queryID     = fastaFile.split("/")[-1].split(".")[0]
    querySeq    = [line for line in open(fastaFile) if not ">" in line ]
    clstrFile   = ''.join([str(work_dir), '/Blastp_cluster/', queryID, '.cluster.clstr'])
    clusterFile = ''.join([str(work_dir), '/Blastp_cluster/', queryID, '.cluster'])
    repID = getClusterRep(queryID, clstrFile)
    reorderFasta(repID, queryID, querySeq, clusterFile)

if __name__ == "__main__": main()
