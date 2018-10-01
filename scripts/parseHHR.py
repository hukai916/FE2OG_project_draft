
class HHSuiteOutput (object):
   class Hit (object):
      def __init__(self,line):
         try:
            self.hitNum = int(line[:3])
            self.name = line[4:35]
            self.prob = float(line[35:40])/100.
            self.numColsAligned = int(line[70:74])
         except ValueError:
            msg = 'could not convert tokens {} to hit'
            raise ValueError(msg.format(line.split()))
   def __init__(self,filename):
      self.parseHHR(filename)
   def parseHHR(self,filename):
      self.hits = []
      readingHits = False
      for line in open(filename):
         tokens = line.split()
         if not tokens:
            readingHits = False
         elif tokens[0] == 'Query':
            self.query = tokens[1]
         elif tokens[0] == 'Match_columns':
            self.numSites = int(tokens[1])
         elif line.startswith(' No Hit'):
            readingHits = True
         elif readingHits:
            self.hits.append(self.Hit(line))

import sys
hhr = HHSuiteOutput(sys.argv[1])
for hit in hhr.hits:
   print hit.name, hit.prob, hit.numColsAligned
