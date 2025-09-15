# 所有数据集
# import sys
# file=sys.argv[1]
# res=sys.argv[2]
# f=open(file)
# res=open(res,'w')
# for i in f:
#     a=i.strip().split('\t')
#     mnvid=a[5]
#     chr=a[0]
#     pos=a[1].split(',')[0]
#     res.write('\t'.join([mnvid,chr,pos])+'\n')
# f.close()
# res.close()


# 单个数据集
import sys
file=sys.argv[1]
dataset=sys.argv[2]
res=sys.argv[3]

f=open(file)
res=open(res,'w')
for i in f:
    a=i.strip().split('\t')
    if dataset in a[-1]:
        mnvid=a[5]
        chr=a[0]
        pos=a[1].split(',')[0]
        res.write('\t'.join([mnvid,chr,pos])+'\n')
f.close()
res.close()