import sys
filename=sys.argv[1] # '/data/jinww/mnv/analyse2/02data/02annotation/gene_canonical.oneline.txt'
res_file=sys.argv[2] # '/data/jinww/mnv/analyse/part4/function/data/'

f=open(filename)
res=open(res_file+'AA_perGene.txt','w')
for i in f:
    a=i.strip().split('\t')
    if "@" in a[-1]:
        b=a[-1].split(' ')
        res.write(a[2]+'\t'+b[1]+'\n')
f.close()
res.close()

# stop_gain
f=open(filename)
res=open(res_file+'stop_change_perGene.txt','w')
for i in f:
    a=i.strip().split('\t')
    if "stop" in a[-1]:
        b=a[-1].split(' ')
        res.write(a[2]+'\t'+b[1]+'\n')
f.close()
res.close()

# stop_lost
f=open(filename)
res=open(res_file+'start_change_perGene.txt','w')
for i in f:
    a=i.strip().split('\t')
    if "start" in a[-1]:
        b=a[-1].split(' ')
        res.write(a[2]+'\t'+b[1]+'\n')
f.close()
res.close()
