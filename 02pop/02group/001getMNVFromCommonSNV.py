import sys

file_common=sys.argv[1] #  file_common='/data/jinww/mnv/analyse2/05population/01overlap/venn/venn_1000G_common.txt'
file_origin=sys.argv[2] # file_origin='/data/jinww/mnv/analyse2/02data/02datasets/1000G.txt'

file2=sys.argv[3] # file2='01origin/mnv.txt'
res=sys.argv[4] # res='02adjust/commonMNV.txt'

commonMNV=[]
with open(file_common) as f:
    for i in f:
        commonMNV.append(i.strip().split('\t')[1])
commonMNV=set(commonMNV)

commonMNV2={}
with open(file_origin) as f:
    for i in f:
        a=i.strip().split('\t')
        if a[2] in commonMNV:
            newa=[]
            chrid=a[0]
            pos=a[1].split(',')
            ref=a[3].split(',')
            alt=a[4].split(',')
            for j in range(len(pos)):
                newa.append(chrid+':'+pos[j]+':'+ref[j]+':'+alt[j])
            newa=','.join(newa)
            commonMNV2[newa]=a[2]

with open(file2) as f, open(res,'w') as w:
    for i in f:
        a=i.strip().split('\t')
        if a[5] in commonMNV2:
            a[2]=commonMNV2[a[5]]
            w.write('\t'.join(a)+'\n')

