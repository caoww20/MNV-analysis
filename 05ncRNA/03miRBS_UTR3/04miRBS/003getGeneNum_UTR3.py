import sys
filename=sys.argv[1] # '/data/jinww/mnv/analyse2/02data/02annotation/gene.oneline.txt'
res_file=sys.argv[2] # 'UTR3_MNVNum.txt'

with open(filename) as f, open(res_file,'w') as res:
    for i in f:
        a=i.strip().split('\t')
        if "UTR3" in a[-1]:
            b=a[-1].split(' ')
            res.write(a[2]+'\t'+b[0]+'\n')
