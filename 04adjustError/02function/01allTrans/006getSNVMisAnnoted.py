#  python 006getSNVMisAnnoted.py merge_total.txt
import sys
file1=sys.argv[1] #file1='merge_total.txt'

snv=[]
diff_snv=[]
with open(file1,'r') as f:
    next(f)
    for line in f:
        line=line.strip().split('\t') 
        mnv_aa=line[14].split('/')[1]

        snv1_id=':'.join(line[15].split(':')[0:-1])
        snv.append(snv1_id)
        snv1_aa=line[18].split('/')[1]
        if snv1_aa != mnv_aa:
            diff_snv.append(snv1_id)
        snv2_id=':'.join(line[19].split(':')[0:-1])
        snv.append(snv2_id)
        snv2_aa=line[22].split('/')[1]
        if snv2_aa != mnv_aa:
            diff_snv.append(snv2_id)
        if line[23] == '.':
            continue
        snv3_id=':'.join(line[23].split(':')[0:-1])
        snv.append(snv3_id)
        snv3_aa=line[26].split('/')[1]
        if snv3_aa != mnv_aa:
            diff_snv.append(snv3_id)
snv=list(set(snv))
diff_snv=list(set(diff_snv))
print(len(snv))
print(len(diff_snv))
print(len(diff_snv)/len(snv))
        