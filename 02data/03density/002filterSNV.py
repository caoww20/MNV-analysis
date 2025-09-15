# python 02filterSNV.py snv_1000G.txt hg19 hg19-blacklist.v2.bed snv_1000G_filter.txt
import sys
file=sys.argv[1] # file='/data/jinww/mnv/analyse/part2/data/density/snv_UKB50w.txt'
MHC=sys.argv[2] # MHC='hg19'
black=sys.argv[3] # black='/data/jinww/mnv/analyse/part2/data/density/hg38-blacklist.v2.bed'
out_file=sys.argv[4] # out_file='/data/jinww/mnv/analyse/part2/data/density/snv_UKB50w_filter.txt'
# hg19 chr6:28477797-33448354
# hg38 chr6:28510120-33480577

black_region={}
for x in range(1,23):
    black_region[str(x)]=[]

f=open(black)
for i in f:
    a=i.strip().split('\t')
    if a[0]!='chrM' and a[0]!='chrY' and a[0]!='chrX':
        black_region[a[0].replace('chr','')].append(a[1]+'-'+a[2])
f.close()

if MHC=='hg19':
    black_region['6'].append('28477797-33448354')
else:
    black_region['6'].append('28510120-33480577')

output=open(out_file,'w')

f=open(file)
for i in f:
    a=i.strip().split('\t')
    flag=0
    for j in black_region[a[1]]:
        if int(a[2])>=int(j.split('-')[0]) and int(a[2])<=int(j.split('-')[1]):
            flag=1
            break
    if flag==0:
        output.write(i) 
f.close()
output.close()