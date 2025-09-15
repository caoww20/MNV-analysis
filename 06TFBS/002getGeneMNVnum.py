# python 002getGeneMNVnum.py protein.promotor.oneline.fix geneMNVNum.txt
import sys
input_file = sys.argv[1] # protein.promotor.oneline.fix
output_file = sys.argv[2]

gene={}
with open(input_file) as f:
    for line in f:
        line = line.strip().split('\t')
        info=line[-1].split(' ')
        geneid=info[6]
        mnvid=line[2]
        if geneid in gene:
            gene[geneid].append(mnvid)
        else:
            gene[geneid]=[mnvid]

for i in gene:
    gene[i]=list(set(gene[i]))

with open(output_file,'w') as w:
    for i in gene:
        w.write(i+'\t'+str(len(gene[i]))+'\n')