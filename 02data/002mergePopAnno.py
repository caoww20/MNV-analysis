# python 03mergePopAnno.py 02adjustOrigin/hg38_humanMNV 03annotation/gene.txt 0 02adjustOrigin/hg38_humanMNV_addgene
import sys
pop_file=sys.argv[1] # 02adjustOrigin/hg38_humanMNV
anno_file=sys.argv[2] # 03annotation/gene.txt
flag=sys.argv[3] # 1 means oneline  
res=sys.argv[4] # 02adjustOrigin/hg38_humanMNV_addgene

# 读入文件
mnv_pop={}
with open(pop_file) as f:
    for line in f:
        line=line.strip().split('\t')
        mnv_pop[line[2]]=line[-1]

# 读入注释文件
with open(anno_file) as f, open(res,'w') as w:
    for line in f:
        line=line.strip().split('\t')
        if flag=='0':
            w.write('\t'.join(line[0:6])+'\t'+mnv_pop[line[2]]+'\t'+line[-1]+'\n')
        else:
            w.write('\t'.join(line[0:6])+'\t'+mnv_pop[line[2]]+'\t'+line[-1].replace(' ','\t')+'\n')