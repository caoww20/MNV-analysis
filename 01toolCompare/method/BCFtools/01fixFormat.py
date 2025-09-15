# 读入参数
import sys
ref_url=sys.argv[1]
file_url=sys.argv[2]
output=sys.argv[3]

pos_dict={}
f=open(ref_url, 'r')
for i in f:
    if 'BCSQ' in i and '#' not in i:
        a=i.strip().split('\t')
        pos_dict[a[1]]=a[3]+','+a[4]
f.close()

#chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
# 打开输入文件a.txt和输出文件b.txt
with open(file_url, 'r') as input_file, open(output, 'w') as output_file:
    for line in input_file:
        a=line.strip().split('\t')
        chrid=a[0]
        pos1=a[1]
        pos1_codon=a[3]+','+a[4]
        pos2=a[7].split('@')[1].split(',')[0]
        pos2_codon=pos_dict[pos2]
        if int(pos1)<int(pos2):
            pos=pos1+','+pos2
            ref=pos1_codon.split(',')[0]+','+pos2_codon.split(',')[0]
            alt=pos1_codon.split(',')[1]+','+pos2_codon.split(',')[1]
        else:
            pos=pos2+','+pos1
            ref=pos2_codon.split(',')[0]+','+pos1_codon.split(',')[0]
            alt=pos2_codon.split(',')[1]+','+pos1_codon.split(',')[1]    

        b=[chrid,pos,'.',ref,alt,'.',str(abs(int(pos1)-int(pos2))),'2','.','.','.']
        output_file.write('\t'.join(b)+'\n')

