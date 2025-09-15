# python 010getInfo.py    
import sys

sig=sys.argv[1] # sig='03fst/filter/sigAFR.txt'
vcf=sys.argv[2] #  vcf='02adjust/MNVFromCommonSNVFilter.vcf'
anno=sys.argv[3] # file_origin='/data/jinww/mnv/analyse2/02data/02annotation/gene.txt'
res=sys.argv[4] # res='03fst/significant/sigAFR.gene'

id={}
for i in range(1,23):
    id[str(i)]=[]

with open(sig) as f:
    next(f)
    for i in f:
        a=i.strip().split('\t')
        id[a[0]].append(a[1])

mnvid={}
for i in range(1,23):
    mnvid[str(i)]=[]  
with open(vcf) as f:
    for i in f:
        if i.startswith('#'):
            continue
        a=i[0:1000].split('\t')
        if a[1] in id[a[0]]:
            mnvid[a[0]].append(a[2])
        

with open(anno) as f, open(res,'w') as w:
    for i in f:
        a=i[0:200].split('\t')
        if a[2] in mnvid[a[0]]:
            w.write(i)