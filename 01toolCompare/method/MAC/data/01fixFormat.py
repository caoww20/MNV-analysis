# 读入参数
import sys
file_url=sys.argv[1]
output=sys.argv[2]

#chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
# 打开输入文件a.txt和输出文件b.txt
with open(file_url, 'r') as input_file, open(output, 'a+') as output_file:
    for line in input_file:
        if 'Mutant' in line:
            a=line.strip().split('\t')
            joint_status=set(a[2])
            if '-' not in joint_status and '0' not in joint_status:
                b=a[0].split(',')
                pos=[]
                ref=[]
                alt=[]
                for i in range(len(b)):
                    #chr22.10519071.G.A
                    pos.append(b[i].split('.')[1])
                    ref.append(b[i].split('.')[2])
                    alt.append(b[i].split('.')[3])
                if abs(int(pos[0])-int(pos[-1]))<=10:
                    b=[file_url.split('_')[0],'22',','.join(pos),'.',','.join(ref),','.join(alt),'.',str(abs(int(pos[0])-int(pos[-1]))),str(len(b)),'.','.','.']
                    output_file.write('\t'.join(b)+'\n')

