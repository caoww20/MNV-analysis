# python 004alluvial.py UTR3/diff_MNV_vs_SNV.txt res/alluvial.txt
import sys
input1=sys.argv[1] # input1='lncRNA/diff_MNV_vs_SNV.txt'
res_file=sys.argv[2] # res_file='res/lncRNA_alluvial.txt'

# 对结果进行去重
mydata=[]
with open(input1) as f:
    for i in f:
        a=i.strip().split('\t')
        id='|'.join(a[0].split('|')[:-1])
        b=id+'\t'+'\t'.join(a[1:])
        mydata.append(b)
mydata=list(set(mydata))

# 考虑no_change的情况
with open(res_file,'w') as res:
    for i in mydata:
        a=i.strip().split('\t')
        mnv_type=a[2]
        if a[-1]=='no_change':
            mnv_type='no_change'
        snv_type=a[3].split(',')
        for ii in snv_type:
            res.write(ii+'\t'+mnv_type+'\n')

# 不考虑no_change的情况
# with open(input1) as f, open(res_file,'w') as res:
#     for i in f:
#         a=i.strip().split('\t')
#         mnv_type=a[2]
#         snv_type=a[3].split(',')
#         for ii in snv_type:
#             res.write(ii+'\t'+mnv_type+'\n')
        