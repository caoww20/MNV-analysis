# python 004getDiff_ms.py snv_diff.res mnv_diff.res merge_res.txt
import sys
file1=sys.argv[1]
file2=sys.argv[2]
output=sys.argv[3]

ref={}
with open(file1) as f1:
    for i in f1:
        i=i.strip().split('\t')
        ref[i[0]]=i[1:]

with open(file2) as f2, open(output,'w') as w:
    for i in f2:
        a=i.strip().split('\t')
        id=a[0].split('|')
        snvid=id[0].split(',')
        snv_GC=[]
        snv_p=[]
        snv_flag=[]
        for ii in snvid:
            snv_GC.append(ref['|'.join([ii]+id[1:])][0])
            snv_p.append(ref['|'.join([ii]+id[1:])][1])
            snv_flag.append(ref['|'.join([ii]+id[1:])][2])
        snv_GC=','.join(snv_GC)
        snv_p=','.join(snv_p)
        snv_flag=','.join(snv_flag)
        change_flag='0'
        if a[-1]=='1' and '1' not in snv_flag:
            change_flag='1'
        newa=a+[snv_GC,snv_p,snv_flag,change_flag]
        w.write('\t'.join(newa)+'\n')

