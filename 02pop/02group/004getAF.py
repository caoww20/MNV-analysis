import sys

AF_file=sys.argv[1] #  file_common='/data/jinww/mnv/analyse2/02data/02adjustOrigin/hg38_humanMNV'
file=sys.argv[2] # file_origin='02adjust/MNVFromCommonSNV.vcf'
res=sys.argv[3] # res='02adjust/MNVFromCommonSNV.freq'

mnvid=[]
with open(file,'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            mnvid.append(line[0:200].split('\t')[2])


with open(AF_file,'r') as f, open(res,'w') as f2:
    for line in f:
        a=line.strip().split('\t')
        if a[2] in mnvid:
            b=a[-1].split('|')
            for i in b:
                if i.startswith('1000G'):
                    c=i.split(';')[1:] # ['AFR:66/3.70e-02','AMR:36/3.67e-02']
                    for j in c:
                        f2.write(a[2]+'\t'+j.replace(':','\t').replace('/','\t')+'\n')
