# python 003getDiff_ms.py 02RNAFold/snv/snv_diff.res 02RNAFold/mnv/mnv_diff.res  02RNAFold/merge_res.txt
import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
output = sys.argv[3]
# 判断0.1 和[0.2,0.3,0.5]的差异，当0.1与所有元素的差值的绝对值都大于1时，返回1，否则返回0
def judge(a,b):
    for i in b:
        if abs(float(a)-float(i))<1:
            return '0'
    return '1'

ref={}
with open(input1) as f1:
    for i in f1:
        i=i.strip().split('\t')
        ref[i[0]]=i[1:]

with open(input2) as f2, open(output,'w') as w:
    for i in f2:
        a=i.strip().split('\t')
        id=a[0].split('|')
        snvid=id[0].split(',')
        snv_MFE=[]
        for ii in snvid:
            snv_MFE.append(ref['|'.join([ii]+id[1:])][1])
        flag=judge(a[2],snv_MFE)
        snv_MFE=','.join(snv_MFE)
        newa=a[0:3]+[snv_MFE,flag]+a[3:]
        w.write('\t'.join(newa)+'\n')

