VioC_18 = [123,125,144,147,
           149,160,161,162,
           173,288,290,295,
           296,297,309,311,
           313,315]

pdbList = ["VioC_VO.pdb", "ANS_2brt.pdb", "AsqJ_5dax.pdb", "CarC_4oj8.pdb",
          "EFE_5lun.pdb", "H6H.pdb", "NapI.pdb", "SadA_3w21.pdb", "SnoN_5equ.pdb",
          "T7H_5c3p.pdb", "WelO5_5iqt.pdb", "AlKB_4qkd.pdb", "CAS_1ds1.pdb",
          "CytC3_3gjb.pdb", "EctD_4q5o.pdb", "LolO.model.pdb", "OrfP_4m25.pdb",
          "SnoK_5epa.pdb", "SyrB2_2fct.pdb", "tauD_1os7.pdb"]

import os
homedir = os.path.dirname(os.getcwd())
pdbFile = os.path.join(homedir, 'PDB/pdbList.txt')
pdbList = [line.split()[0] for line in open(pdbFile)]
tmScoreFile = os.path.join(homedir, 'PDB/tm_score.txt')


def filterG(g): # return g with limited number of nodes
    for key in g:
        hitList = g[key]
        tem = []
        if len(hitList) > 5:
            scoreList = [x.split("+")[1] for x in hitList]
            candidates = sorted(map(float, scoreList))[:]
            for item in hitList:
                if float(item.split("+")[1]) in candidates:
                    tem.append(item.split("+")[0])
        else:
            for item in hitList:
                tem.append(item.split("+")[0])
        g[key] = tem
    return(g)

g = {}
for line in open(tmScoreFile):
    tem = line.split()
    if not tem[0] in g:
        g[tem[0]] = []
        if float(tem[2]) > 0.5:
            g[tem[0]].append(tem[1] + "+" + tem[2])
    else:
        if float(tem[2]) > 0.5:
            g[tem[0]].append(tem[1] + "+" + tem[2])

g = filterG(g)


if __name__ == "__main__":
    print(pdbList)
    g = filterG(g)
    for key in g:
        print(key, g[key])
"""
'AsqJ_5dax.pdb', 'OrfP_4m25.pdb', 'CarC_4oj8.pdb', 'SyrB2_2fct.pdb',
"""
