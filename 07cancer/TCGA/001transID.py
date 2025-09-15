# python 001transID.py TCGA.mnv gene.txt eQTL.txt eQTL_fix.txt
import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]
output = sys.argv[4]

mnvid={}
with open(input1) as f:
    for line in f:
        line = line.strip().split('\t')
        id='_'.join(line[0:2]+line[3:6])
        info=line[-1]
        af=''
        for i in info.split('|'):
            if i.startswith('TCGA'):
                af=i
        mnvid[id]=[line[2],af]

gene={}
with open(input2) as f:
    for line in f:
        line = line.strip().split('\t')
        gene[line[2]]=line[-1]

with open(input3) as f,open(output,'w') as w:
    for line in f:
        i=line
        line = line.strip().split('\t')
        if line[2] in mnvid:
            info=mnvid[line[2]][1]
            line[2]=mnvid[line[2]][0]
            w.write('\t'.join(line+[gene[line[2]],info])+'\n')