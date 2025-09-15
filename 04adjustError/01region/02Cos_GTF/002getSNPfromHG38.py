import sys
input=sys.argv[1]
type=sys.argv[2] # gene or other or total
flag=sys.argv[3] # utr3 or utr5 or exon or splice or intron
output=sys.argv[4]
        
# 读入humanMNV的数据库
with open(input) as f, open(output,'w') as res:
    for i in f:
        a=i.strip().split('\t')
        chr=a[0]
        pos_list=a[1].split(',')
        if a[-1] !='.':
            if type=='gene':
                if flag in a[-1]:
                    for ii in pos_list:
                        res.write(chr+'\t'+ii+'\n')
            else:
                for ii in pos_list:
                    res.write(chr+'\t'+ii+'\n')