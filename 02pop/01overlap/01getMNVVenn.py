import sys
file=sys.argv[1] # /data/jinww/mnv/analyse2/02data/02adjustOrigin/hg38_humanMNV
resname=sys.argv[2] # ./venn/venn_1000G.txt

# pop
mydata2={}
mydata2['AFR']=[]
mydata2['AMR']=[]
mydata2['EAS']=[]
mydata2['EUR']=[]
mydata2['SAS']=[]

# def getMAFMNV(x,cancer):
#     for ii in x.split(';'):
#         if ii.startswith(cancer):
#             maf=float(ii.split('/')[1].split('|')[0])
#             if maf>0.5:
#                 maf=1-maf
#             return maf


# 1000G
f=open(file)
for i in f:
    a=i.strip().split('\t')
    mnvid=a[2]
    # Count presence per population. Exclude combined tags like EUR-AFR, hence we check for ';POP:'
    if ';AFR:' in a[-1]:
        mydata2['AFR'].append(mnvid)
    if ';AMR:' in a[-1]:
        mydata2['AMR'].append(mnvid)
    if ';EAS:' in a[-1]:
        mydata2['EAS'].append(mnvid)
    if ';EUR:' in a[-1]:
        mydata2['EUR'].append(mnvid)
    if ';SAS:' in a[-1]:
        mydata2['SAS'].append(mnvid)

f.close()


res=open(resname,'w')
for i in mydata2:
    for ii in mydata2[i]:
        res.write(i+'\t'+ii+'\n')
res.close()

