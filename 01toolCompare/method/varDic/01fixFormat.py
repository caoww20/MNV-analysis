# 读入参数
import sys
file_url=sys.argv[1]
output=sys.argv[2]

#chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
# 打开输入文件a.txt和输出文件b.txt
with open(file_url, 'r') as input_file, open(output, 'w') as output_file:
    for line in input_file:
        if '#' not in line:  # 如果行不包含 #
            a=line.strip().split('\t')
            pos1=int(a[1])
            ref=list(a[3])
            alt=list(a[4])
            if len(ref)==len(alt) and len(ref)<=11:
                pos=[]
                for i in range(len(ref)):
                    pos.append(str(pos1+i))
                b=[file_url.split('/')[0],'22',','.join(pos),'.',','.join(ref),','.join(alt),'.',str(len(ref)-1),str(len(ref)),'.','.','.']
                output_file.write('\t'.join(b)+'\n')