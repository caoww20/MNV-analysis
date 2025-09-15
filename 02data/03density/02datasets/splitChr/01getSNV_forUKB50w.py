## 获得SNV信息
import sys,gzip
file=sys.argv[1]
res=sys.argv[2]
def openfile(filename, mode="r"):
    if filename[-3:] == '.gz':
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

res=open(res,'w')
with openfile(file) as f:
    for i in f:
        i = i[0:3000].decode('utf-8')
        if i[0]!='#':
            a=i[0:2000].split('\t')
            chr=a[0]
            pos=a[1]
            snv=a[2]
            if len(a[3])==1 and len(a[4])==1:
                res.write('\t'.join([snv,chr,pos])+'\n')
res.close()