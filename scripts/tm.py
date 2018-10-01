from tmalign import tm_align_all, parse_tm_output
from Constant import pdbList

res = tm_align_all(pdbList)
for x in res:
    print("\t".join(map(str, [x[0], x[1], x[2]])))
