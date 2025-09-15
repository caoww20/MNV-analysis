# python 003getDiff_mnv.py 02RNAFold/ref/ref_fix.res 02RNAFold/mnv/mnv_fix.res 02RNAFold/mnv/mnv_diff.res 
import sys
input1 = sys.argv[1]
input2 = sys.argv[2]
output = sys.argv[3]

ref={}
with open(input1) as f1:
    for i in f1:
        i=i.strip().split('\t')
        id=i[0].split('|')
        snvid=id[0].split(',')
        for ii in snvid:
            ref['|'.join([ii]+id[1:])]=i[1:]

with open(input2) as f2, open(output,'w') as w:
    for i in f2:
        # id MFE diff_MFE ref_seq ref_struc mnv_seq mnv_struc
        a=i.strip().split('\t')
        diff_MFE=str(round(float(a[1])-float(ref[a[0]][0]),2))
        newa=[a[0],a[1],diff_MFE,ref[a[0]][1],ref[a[0]][2],a[2],a[3]]
        w.write('\t'.join(newa)+'\n')

