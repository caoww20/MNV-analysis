import os,sys
#  python 005mergeFST.py ./03fst/ AFR,AMR,EAS,EUR,SAS 0 single_merge & 

# path = "/data/jinww/mnv/analyse2/05population/02group/03fst/" # 路径
path=sys.argv[1]
file_list = os.listdir(path) # 获取路径下所有文件名
file_list = [f for f in file_list if "log" not in f and "txt" not in f] # 过滤文件名包含"fst"的文件
pop =sys.argv[2]
pop=pop.split(',')
flag=sys.argv[3] # 0表示single 1表示region
res=sys.argv[4]

# 第一次选择文件
fst_select1=[]
for i in file_list:
    for ii in pop:
        if ii in i:
            fst_select1.append(i)
            break
# 第二次选择文件
fst_files=[]
for i in fst_select1:
    if flag=='0':
        if '100k' not in i:
            fst_files.append(i)
    else:
        if '100k' in i:
            fst_files.append(i)

# 读取文件
if flag=='0':
    info={}
    res=res+'.txt'
    with open(path+res,'w') as w:
        header=['chrom','pos']+[i.split('.')[0] for i in fst_files]
        w.write('\t'.join(header)+'\n')

    with open(path+fst_files[0], "r") as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = line[1]
            fst = line[2]
            info[chrom+"_"+pos]=[fst]

    for fst_file in fst_files[1:]:
        with open(path+fst_file, "r") as f:
            next(f)
            for line in f:
                line = line.strip().split('\t')
                chrom = line[0]
                pos = line[1]
                fst = line[2]
                if chrom+"_"+pos in info:
                    info[chrom+"_"+pos].append(fst)
    with open(path+res,'a+') as w:
        for k,v in info.items():
            v2='\t'.join(v)
            if "-" not in v2 and "nan" not in v2:
                w.write(k.replace('_','\t')+'\t'+v2+'\n')
else:
    res1=res+'_weight.txt'
    res2=res+'_mean.txt'

    with open(path+res1,'w') as w:
        header=['chrom','start','end']+[i.split('.')[0] for i in fst_files]
        w.write('\t'.join(header)+'\n')
    with open(path+res2,'w') as w:
        header=['chrom','start','end']+[i.split('.')[0] for i in fst_files]
        w.write('\t'.join(header)+'\n')

    info={}
    with open(path+fst_files[0], "r") as f:
        next(f)
        for line in f:
            line = line.strip().split('\t')
            chrom = line[0]
            start = line[1]
            end = line[2]
            fst = ','.join(line[4:])
            info[chrom+"_"+start+"_"+end]=[fst]

    for fst_file in fst_files[1:]:
        with open(path+fst_file, "r") as f:
            next(f)
            for line in f:
                line = line.strip().split('\t')
                chrom = line[0]
                start = line[1]
                end = line[2]
                fst = ','.join(line[4:])
                if chrom+"_"+start+"_"+end in info:
                    info[chrom+"_"+start+"_"+end].append(fst)
    
    with open(path+res1,'a+') as w:
        for k,v in info.items():
            v2='\t'.join([i.split(',')[0] for i in v])
            if "-" not in v2 and "nan" not in v2:
                w.write(k.replace('_','\t')+'\t'+v2+'\n')
    with open(path+res2,'a+') as w:
        for k,v in info.items():
            v2='\t'.join([i.split(',')[1] for i in v])
            if "-" not in v2 and "nan" not in v2:
                w.write(k.replace('_','\t')+'\t'+v2+'\n')

