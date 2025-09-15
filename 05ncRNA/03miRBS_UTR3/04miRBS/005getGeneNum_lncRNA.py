import sys
filename=sys.argv[1] # '/data/jinww/mnv/analyse3/02data/02annotation/lnc.oneline.txt'
res_file=sys.argv[2] # 'lncRNA_MNVNum.txt'

with open(filename) as f, open(res_file,'w') as res:
    for i in f:
        a=i.strip().split('\t')
        if "all" in a[-1] and 'exon' in a[-1] and 'intro' not in a[-1]:
            b=a[-1].split(' ')
            res.write(a[2]+'\t'+b[8]+'\n')
