# python 007getFSTSH.python MNVFromCommonSNV.vcf AFR,AMR,EAS,EUR,SAS a.sh
import sys
file=sys.argv[1] # MNVFromCommonSNV.vcf
pop=sys.argv[2]
res=sys.argv[3]
pop=pop.split(',')
# 将pop两两组合形成新的列表
pop2=[]
for i in range(len(pop)):
    for j in range(i+1,len(pop)):
        pop2.append(pop[i]+'_'+pop[j])
# 输出结果
with open(res,'w') as w:
    a=[]
    for i in range(len(pop)):
        a.append(" --weir-fst-pop 01origin/sample/"+pop[i]+".sample")
    a=''.join(a)
    w.write("vcftools --vcf 02adjust/"+file+a+" --out 03fst/all &\n")
    w.write("vcftools --vcf 02adjust/"+file+a+" --fst-window-size 100000 --fst-window-step 10000 --out 03fst/all_100k &\n")

with open(res,'a+') as w:
    for i in pop2:
        a="vcftools --vcf 02adjust/"+file+" --weir-fst-pop 01origin/sample/"+i.split('_')[0]+".sample --weir-fst-pop 01origin/sample/"+i.split('_')[1]+".sample --out 03fst/"+i+" &"
        w.write(a+'\n')
        b="vcftools --vcf 02adjust/"+file+" --weir-fst-pop 01origin/sample/"+i.split('_')[0]+".sample --weir-fst-pop 01origin/sample/"+i.split('_')[1]+".sample --fst-window-size 100000 --fst-window-step 10000 --out 03fst/"+i+"_100k &"
        w.write(b+'\n')
