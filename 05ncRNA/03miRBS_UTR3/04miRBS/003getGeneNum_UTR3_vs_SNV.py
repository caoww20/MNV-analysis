# python 003getGeneNum_UTR3_vs_SNV.py /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt UTR3/diff_MNV_vs_SNV.txt res/UTR3_MNVNum_vs_SNV.txt
import sys
input1=sys.argv[1] # '/data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt'
input2=sys.argv[2] # 'UTR3/diff_MNV_vs_ref.txt'
res_file=sys.argv[3] # 'UTR3_MNVNum_diff.txt'

gene={}
with open(input1) as f:
    for i in f:
        a=i.strip().split('\t')
        gene[a[0]]=a[3]

with open(input2) as f, open(res_file,'w') as res:
    for i in f:
        a=i.strip().split('\t')
        if a[-1]=='no_change':
            continue
        a=a[0].split('|')[1:]
        mnvid=a[0]
        trans_id=a[1]
        geneid=gene[trans_id]
        res.write(mnvid+'\t'+geneid+'\n')
        
