import sys
file=sys.argv[1] # file='/data/jinww/mnv/analyse/part2/data/snv_GTEx.txt'
flag=sys.argv[2] # hg19 or hg38
out_file=sys.argv[3] #'/data/jinww/mnv/analyse/part2/data/snv_GTEx_Density.txt'
# hg19 chr6:28477797-33448354
# hg38 chr6:28510120-33480577

def getDes(lst,step):
    # lst = set([0, 1, 2, 8, 9, 17, 18])  # 给定的列表

    # step = 3  # 步长
    result = []  # 存储每个5步长内的存在的值的数量

    for i in range(0, max(lst)+1, step):
        count = 0  # 记录当前5步长内存在的值的数量
        for j in lst:
            if j in range(i, i+step):
                count += 1
        result.append(count)

    # 如果步长内没有值则输出0
    for i in range(len(result)):
        if result[i] == 0:
            result[i] = 0

    # print(result)  # 输出新的列表
    return(result)

output=open(out_file,'w')
for chrid in [str(x) for x in range(1,23)]:
    pos=[]
    f=open(file)
    for i in f:
        a=i.strip().split('\t')
        if a[1]==str(chrid):
            if a[1]=='6':
                if flag=='hg19':
                    if int(a[2])<28477797 or int(a[2])>33448354:
                        pos.append(int(a[2]))
                elif flag=='hg38':
                    if int(a[2])<28510120 or int(a[2])>33480577:
                        pos.append(int(a[2]))
            else:
                pos.append(int(a[2]))
    f.close()
    res=getDes(set(pos),1000000)
    for i in range(len(res)):
        output.write(str(chrid)+'\t'+str(i)+'\t'+str(res[i])+'\n')

output.close()