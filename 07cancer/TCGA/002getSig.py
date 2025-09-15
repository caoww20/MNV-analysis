# python 002getSig.py eQTL_fix.txt UKB50wMNVGWAS/03.Result/01.allCancer/14.lung/fastGWA_GLMM_final.fastGWA eQTL_gwas.txt
import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
output = sys.argv[3]


gwas_id={}
with open(input2) as f:
    next(f)
    for line in f:
        line = line.strip().split('\t')
        # beta pvalue
        gwas_id[line[1]]=[line[10],line[12]]

with open(input1) as f, open(output,'w') as res:
    for line in f:
        line = line.strip().split('\t')
        if line[0] in ['LUAD','LUSC']:
            # 如果 ['promotor','UTR5','CDS','splice','UTR3'] 中任一一个字符串在一个长的字符串a中则输出
            if any(x in line[-2] for x in ['up:','UTR5','exon','splice','UTR3']):
                mnvid=line[2]
                if mnvid in gwas_id:
                    res.write('\t'.join(line[:-2])+'\t'+'\t'.join(gwas_id[mnvid]+line[-2:])+'\n')
                else:
                    res.write('\t'.join(line[:-2])+'\t'+'\t'.join(['NA','NA']+line[-2:])+'\n')