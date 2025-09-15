import sys
input=sys.argv[1]
output=sys.argv[2]
        
# 读入humanMNV的数据库
with open(input) as f, open(output,'w') as res:
    for i in f:
        a=i.strip().split('\t')
        chr=a[0]
        pos_list=a[1].split(',')
        for ii in pos_list:
            res.write(chr+'\t'+ii+'\n')