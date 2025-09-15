import sys
filename=sys.argv[1]
resname=sys.argv[2]

## 提取hg38的数据中所有落在同一个密码子上的MNV（2joint MNV distance<=2 or 3joint MNV distance ==2）##########################################
f=open(filename)
res=open(resname,'w')
for i in f:
    a = i.strip('\n').split('\t')
    mnvDistance=int(a[1].split(',')[-1])-int(a[1].split(',')[0])
    ifcodon=len(a[-3].split(','))
    flag=0
    for ii in a[-5].split(','):
        if 'exon' not in ii:
            flag=1
    if mnvDistance<=2 and ifcodon==1 and flag==0: # 确保在一个密码子内
        res.write(i)      
f.close()
res.close()