# python 003getGene.py eQTL_fix.txt eQTL_fix.gene
import sys,re
input1 = sys.argv[1]
output = sys.argv[2]

# Split by mutiple symbols
def go_split(s, symbol):
    symbol = "[" + symbol + "]+"
    result = re.split(symbol, s)
    return [x for x in result if x]

with open(input1) as f, open(output,'w') as res:
    gene=[]
    for line in f:
        line = line.strip().split('\t')
        info=line[-2]
        if info!='.':
            tmp=go_split(info,'$|')
            for ii in tmp:
                gene.append(ii.split(' ')[1])
    gene=set(gene)
    for i in gene:
        res.write(i+'\n')