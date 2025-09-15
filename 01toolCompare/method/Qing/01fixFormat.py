# 读入参数
import sys
file_url=sys.argv[1]
output=sys.argv[2]

#chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
# 打开输入文件a.txt和输出文件b.txt
with open(file_url, 'r') as input_file, open(output, 'w') as output_file:
    for line in input_file:
        if 'locus' not in line:  # 如果行不包含 #
            a=line.strip().split('\t')
            chrid=a[0].split(':')[0]
            pos1=a[0].split(':')[1]
            pos2=a[2].split(':')[1]
            pos1_codon=a[1].replace('[','').replace(']','').replace('\"','')
            pos2_codon=a[3].replace('[','').replace(']','').replace('\"','')
            if int(pos1)<int(pos2):
                pos=pos1+','+pos2
                ref=pos1_codon.split(',')[0]+','+pos2_codon.split(',')[0]
                alt=pos1_codon.split(',')[1]+','+pos2_codon.split(',')[1]
            else:
                pos=pos2+','+pos1
                ref=pos2_codon.split(',')[0]+','+pos1_codon.split(',')[0]
                alt=pos2_codon.split(',')[1]+','+pos1_codon.split(',')[1]        
            b=[chrid,pos,'.',ref,alt,'.',str(abs(int(pos1)-int(pos2))),'2',a[6],'.',a[5]]
            output_file.write('\t'.join(b)+'\n')