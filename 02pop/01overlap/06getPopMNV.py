import sys
file=sys.argv[1]
pop=sys.argv[2]
res=sys.argv[3]

# pop id chr pos
res=open(res,'w')
f=open(file)
for i in f:
    a=i.strip().split('\t')
    mnvid=a[5]
    chr=a[0]
    pos=a[1].split(',')[0]
    res.write('\t'.join([pop,mnvid,chr,pos])+'\n')
f.close()
res.close()
