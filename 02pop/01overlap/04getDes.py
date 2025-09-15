import sys
file=sys.argv[1] # file='/data/jinww/mnv/analyse/part2/data/snv_GTEx.txt'
out_file=sys.argv[2] #'/data/jinww/mnv/analyse/part2/data/snv_GTEx_Density.txt'

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
        if a[2]==str(chrid):
            pos.append(int(a[3]))
    f.close()
    res=getDes(set(pos),1000000)
    for i in range(len(res)):
        output.write(str(chrid)+'\t'+str(i)+'\t'+str(res[i])+'\n')

output.close()