import math,sys
filename=sys.argv[1] # filename='gnomAD_fix_snv_anno.oneline.txt'
anno=sys.argv[2] # anno='/data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt'
resname=sys.argv[3] # resname='gnomAD_fix.snv.new'

with open(filename) as f, open(resname,'w') as w:
    for line in f:
        a = line.strip('\n').replace(' ','\t').split('\t')
        newa= [a[0],a[1],a[3],a[4]]+a[6:-1]
        if 'exon' in a[9]:
            w.write('\t'.join(newa)+'\n')

## 提取snv的3碱基##########################################
# 根据经典转录本将其CDS序列从en_GRCh38.108_canonical_gene.anno.txt提取出来
f=open(anno)
mytrans={}
for i in f:
    a = i.strip('\n').split('\t')
    mytrans[a[0]]=a[20]
f.close()
# 读入snv的注释数据并提取snv结果
f=open(resname)
s=f.readlines()
f.close()
f1=open(resname,'w')
for i in s:
    a = i.strip('\n').split('\t')
    pos=int(a[8][1:-1])
    aa_pos = math.ceil(pos/3)
    ref=mytrans[a[6]][((aa_pos-1)*3):aa_pos*3].lower()
    c1=ref[0]
    c2=ref[1]
    c3=ref[2]
    flag=pos-1-(aa_pos-1)*3
    if flag==0:
        alt=a[8][-1]+c2+c3
    elif flag==1:
        alt=c1+a[8][-1]+c3
    else:
        alt=c1+c2+a[8][-1]
    newa=a+[ref+'/'+alt,a[9][0]+'/'+a[9][-1]]
    f1.write('\t'.join(newa)+'\n')
f.close()
f1.close()