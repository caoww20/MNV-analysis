# python 003getSig.py res_RNAsnp.txt ./
import sys
file=sys.argv[1]
url=sys.argv[2]

with open(file, 'r') as f, open(url+'mnv_diff.res', 'w') as w1, open(url+'snv_diff.res', 'w') as w2:
    for i in f:
        line=i.strip().split('\t')
        id=line[1]
        CG=line[4]
        window=int(line[2])
        if window<500: # 等价于 if RNAsnp_mode==1:
            p_value=line[7]
        else:
            p_value=line[-1]
        flag='0'
        if float(p_value)<=0.2:
            flag='1'
        newa=[id,CG,p_value,flag]
        if '-' in line[0]:
            w1.write('\t'.join(newa)+'\n')
        else:
            w2.write('\t'.join(newa)+'\n')
        
