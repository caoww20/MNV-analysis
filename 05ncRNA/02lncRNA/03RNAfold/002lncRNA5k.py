# python 002lncRNA5k.py lncRNA_MNV.seq lncRNA_MNV_5k.seq
import sys
input1 = sys.argv[1]
output = sys.argv[2]
with open(input1) as f, open(output,'w') as w:
    for line in f:
        if line.startswith('>'):
           a=line
           b=next(f)
           if len(b)<=5000:
              w.write(a)
              w.write(b)