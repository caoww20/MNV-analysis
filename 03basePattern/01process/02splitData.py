import sys
##从原始表格中提取不同连体的MNV数据#######################################################
filename = sys.argv[1] 
# filename='1000G.txt'
output=filename.split('.')[0]+'_adjust.txt'
# 提取MNV矫正MNV和原始MNV一致的，因为这部分的MNV具有唯一性
f=open(filename)
res=open(output,'w')
# 22      50786565,50786574       .       G,A     C,G     22:50786565:G:C,22:50786574:A:G 9       2       1       1       0.0001561524    0       3119
# mnvid	refs	alts	distance mnvtype	snp1,snp2	ac1,ac2	ac_mnv
# 1:10247	"T,A"	"C,T"	1	"1-10247-T-C","1-10248-A-T"	23,19	4   
for i in f:
    a=i.strip().split()
    ismulti=len(a[2].split('.'))
    if ismulti==1:
        chrid=a[0]
        if int(a[7])==2:
            snv1=':'.join([chrid,a[1].split(',')[0],a[3].split(',')[0],a[4].split(',')[0]])
            snv2=':'.join([chrid,a[1].split(',')[1],a[3].split(',')[1],a[4].split(',')[1]])
            res.write('\t'.join([a[2],a[3],a[4],a[6],a[7],','.join([snv1,snv2]),'0,0',a[8].split('/')[0].split(':')[1]])+'\n')
        elif int(a[7])==3:
            snv1=':'.join([chrid,a[1].split(',')[0],a[3].split(',')[0],a[4].split(',')[0]])
            snv2=':'.join([chrid,a[1].split(',')[1],a[3].split(',')[1],a[4].split(',')[1]])
            snv3=':'.join([chrid,a[1].split(',')[2],a[3].split(',')[2],a[4].split(',')[2]])
            res.write('\t'.join([a[2],a[3],a[4],a[6],a[7],','.join([snv1,snv2,snv3]),'0,0,0',a[8].split('/')[0].split(':')[1]])+'\n') 
        elif int(a[7])==4:
            snv1=':'.join([chrid,a[1].split(',')[0],a[3].split(',')[0],a[4].split(',')[0]])
            snv2=':'.join([chrid,a[1].split(',')[1],a[3].split(',')[1],a[4].split(',')[1]])
            snv3=':'.join([chrid,a[1].split(',')[2],a[3].split(',')[2],a[4].split(',')[2]])
            snv4=':'.join([chrid,a[1].split(',')[3],a[3].split(',')[3],a[4].split(',')[3]])
            res.write('\t'.join([a[2],a[3],a[4],a[6],a[7],','.join([snv1,snv2,snv3,snv4]),'0,0,0,0',a[8].split('/')[0].split(':')[1]])+'\n')    
        elif int(a[7])==5:
            snv1=':'.join([chrid,a[1].split(',')[0],a[3].split(',')[0],a[4].split(',')[0]])
            snv2=':'.join([chrid,a[1].split(',')[1],a[3].split(',')[1],a[4].split(',')[1]])
            snv3=':'.join([chrid,a[1].split(',')[2],a[3].split(',')[2],a[4].split(',')[2]])
            snv4=':'.join([chrid,a[1].split(',')[3],a[3].split(',')[3],a[4].split(',')[3]])
            snv5=':'.join([chrid,a[1].split(',')[4],a[3].split(',')[4],a[4].split(',')[4]])
            res.write('\t'.join([a[2],a[3],a[4],a[6],a[7],','.join([snv1,snv2,snv3,snv4,snv5]),'0,0,0,0,0',a[8].split('/')[0].split(':')[1]])+'\n')          
f.close()
res.close()


##匹配ac1和ac2######################################################
# 读入rs的信息 
f=open(filename.split('.')[0]+'_snv.txt')
snv={}
for i in f:
    # 1:13526:C:T     1       0.000596659
    a=i.strip().split('\t')
    snv[a[0]]=str(a[1])
f.close()
## 将rs的AC信息合并到数据中 
## 2joint
f=open(output)
s=f.readlines()
f.close()
f1=open(output,'w')
f1.write('mnvid\trefs\talts\tdistance\tmnvtype\tsnvid\tsnv\tmnv\trate\n')
for i in s:
    a=i.strip().split('\t')
    b=[]
    for ii in a[5].split(','):
        b.append(snv[ii])
    a[6]=','.join(b)
    a.append(str(int(a[7])/int(b[0])))
    f1.write('\t'.join(a)+'\n')
f1.close()

