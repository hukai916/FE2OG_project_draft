class FatCatOutput():
    class Hit():
        def __init__(self, line):
            try:
                tokens = line.split()
                self.pdb1 = tokens[0]
                self.pdb1length = int(tokens[1])
                self.pdb2 = tokens[2]
                self.pdb2length = int(tokens[3])
                self.score = float(tokens[4])
                self.pvalue = float(tokens[5])
                self.twist = int(tokens[6])
                self.optlen = int(tokens[7])
                self.optrmsd = float(tokens[8])
                self.chainrmsd = float(tokens[9])
                self.alignlen = int(tokens[10])
                self.gap = int(tokens[11])
            except ValueError:
                msg = 'Can not covert to Hit object'
                raise ValueError(msg)
    def __init__(self, filename):
        self.parseFatcat(filename)
    def parseFatcat(self, filename):
        self.hits = []
        for line in open(filename):
            if not line.startswith("#"):
                #print(line)
                self.hits.append(self.Hit(line))

import sys
filename = sys.argv[1]
fatcat_output = FatCatOutput(filename)
for hit in fatcat_output.hits:
    if hit.optlen / hit.pdb1length > 0.75:
        print(hit.pdb2)
