import sys
filename=sys.argv[1]
resname=sys.argv[2]
## 将提取的MNV转化为SNV格式进行注释##########################################
f=open(filename)
snv=[]
for i in f:
    a = i.strip('\n').split('\t')
    chr=a[0]
    pos=a[1].split(',')
    ref=a[3].split(',')
    alt=a[4].split(',')
    for i in range(len(pos)):
        b=[chr,pos[i],'.',ref[i],alt[i],'.']
        snv.append('\t'.join(b)+'\n')     
f.close()
snv=list(set(snv))
res=open(resname,'w')
for i in snv:
    res.write(i)
res.close()
