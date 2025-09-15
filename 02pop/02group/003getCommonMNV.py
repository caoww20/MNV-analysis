import sys

file_common=sys.argv[1] #  file_common='/data/jinww/mnv/analyse2/05population/01overlap/venn/venn_1000G_common.txt'
vcf=sys.argv[2] # file_origin='02adjust/MNVFromCommonSNV.vcf'
res=sys.argv[3] # res='02adjust/commonMNV.txt'

commonMNV={}
with open(file_common) as f:
    for i in f:
        a=i.strip().split('\t')
        if a[1] not in commonMNV:
            commonMNV[a[1]]=1
        else:
            commonMNV[a[1]]=commonMNV[a[1]]+1

commonMNV2=[]
for i in commonMNV:
    if commonMNV[i]==5:
        commonMNV2.append(i)


with open(vcf) as f, open(res,'w') as w:
    for i in f:
        if i.startswith('#'):
            w.write(i)
        else:
            mnvid=i[0:200].split('\t')[2]
            if mnvid in commonMNV2:
                w.write(i)

