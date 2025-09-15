# 读入参数
import sys
file_url=sys.argv[1]
output=sys.argv[2]

# 将列表[1,2,5,6,10,11]进行拆分，以窗口为2为例，拆分为[[1,2],[5,6],[10,11]]
def split_list(list, window):
    return [list[i:i+window] for i in range(0, len(list), window)]

# split_list([1,2,5,6,10,11], 2)
pos_list=[]
pos_dict={}
f=open(file_url, 'r')
for i in f:
    if 'CHROM' not in i:
        a=i.strip().split('\t')
        pos_list.append(int(a[1]))
        pos_dict[int(a[1])]=[a[3],a[4]]
f.close()
pos_list=sorted(pos_list)
res=split_list(pos_list, 2)

#chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
f=open(output, 'a+')
for i in res:
    b=[file_url.split('_')[0],'22',str(i[0])+','+str(i[1]),'.',pos_dict[i[0]][0]+','+pos_dict[i[1]][0],pos_dict[i[0]][1]+','+pos_dict[i[1]][1],'.',str(i[1]-i[0]),'2','.','.','.']
    f.write('\t'.join(b)+'\n')
f.close()

