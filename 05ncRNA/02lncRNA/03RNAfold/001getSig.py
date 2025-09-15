# python 001getSig.py ../02RNAsnp/merge_res.txt ../01data/lncRNA_adjust.txt res_sig.mnv
import sys
file=sys.argv[1] # file='../02RNAsnp/merge_res.txt'
mnv=sys.argv[2] # mnv='../01data/lncRNA_adjust.txt'
output=sys.argv[3] # output='res_sig.mnv'

mnvid=[]
with open(file, 'r') as f:
    for i in f:
        line=i.strip().split('\t')
        id=line[0].split('|')[1:3]
        if line[-1]=='1':
            mnvid.append('|'.join(id))
mnvid=set(mnvid)

with open(mnv, 'r') as f, open(output, 'w') as w:
    for i in f:
        line=i.strip().split('\t')
        if '|'.join([line[2],line[-1].split(' ')[6]]) in mnvid:
            w.write(i)