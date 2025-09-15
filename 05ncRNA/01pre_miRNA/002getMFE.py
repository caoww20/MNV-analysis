# python 002getMFE.py 02RNAFold/ref/ref.res 02RNAFold/ref/ref_fix.res
import sys,re
input=sys.argv[1]
output=sys.argv[2]

with open(input) as f, open(output,'w') as w:
    # 如果是>开头，则输出，如果是>开头的后一行也输出，其他情况不输出
    for line in f:
        if line.startswith('>'):
            # id=line[1:].strip().replace('_','\t')
            id=line[1:].strip().replace('_','|')
            seq=next(f).strip()
            MFE=next(f)
            struc=MFE.split(' ')[0]
            MFE=re.search(r'-?\d+\.\d+', MFE)
            MFE=MFE.group()
            w.write(id+'\t'+MFE+'\t'+seq+'\t'+struc+'\n')
            