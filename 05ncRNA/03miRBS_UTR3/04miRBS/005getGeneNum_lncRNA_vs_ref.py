# python 005getGeneNum_lncRNA_vs_ref.py /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt lncRNA/diff_MNV_vs_ref.txt res/lncRNA_MNVNum_vs_ref.txt
import sys
input1=sys.argv[1] # '/data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt'
input2=sys.argv[2] # 'UTR3/diff_MNV_vs_ref.txt'
res_file=sys.argv[3] # 'UTR3_MNVNum_diff.txt'

gene={}
with open(input1) as f:
    for i in f:
        a=i.strip().split('\t')
        gene[a[5]]=a[7]

with open(input2) as f, open(res_file,'w') as res:
    for i in f:
        a=i.strip().split('\t')
        if a[2]=='-' or a[-2]=='-' or a[-1]=='unclear':
            continue
        a=a[0].split('|')[1:]
        mnvid=a[0]
        trans_id=a[1]
        geneid=gene[trans_id]
        res.write(mnvid+'\t'+geneid+'\n')
        
