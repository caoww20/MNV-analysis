import numpy as np
import sys

fa_url=sys.argv[1]
cutoff=int(sys.argv[2]) #5
filename=sys.argv[3]
resfile=sys.argv[4]
# python 04addSeq.py /data/jinww/mnv/MNVAnno/database/human_build/ensembl/hg38/chr_fix.fa 5 mnv.txt

# 根据染色体读入参考基因组
f=open(fa_url)
s=f.readlines()
chr_seq = {}
for i in range(0,len(s),2):
    chr_seq[s[i].strip('\n').replace('>','')]=[s[i+1].strip('\n'),len(s[i+1].strip('\n'))]
del s
f.close()

# 进行追加数据
f=open(filename)
s=f.readlines()
f.close()
f1=open(resfile,'w')
f1.write('mnvid\trefs\talts\tdistance\tmnvtype\tsnvid\tsnv\tmnv\trate\tref_seq\talt_seq\n')
s=s[1:]

# 添加序列
for i in s:
    a=i.strip().split('\t')
    b=a[5].split(',')
    chrid=b[0].split(':')[0]
    minPos=int(b[0].split(':')[1])-1 # python中对应的碱基位置
    maxPos=int(b[-1].split(':')[1])-1 # python中对应的碱基位置
    pos=int(np.floor(np.median([minPos,maxPos])))
    ## ref
    ref_seq=chr_seq[chrid][0][(pos-cutoff+1):(pos+cutoff+1)]
    a.append(ref_seq)
    ## alt
    ## 先从seq中提取一段，然后再改变，再拼接
    pos_list=list(range(minPos-20,maxPos+20))
    pos_seq=list(chr_seq[chrid][0][(minPos-20):(maxPos+20)])
    for ii in b:
        fixPos=int(ii.split(':')[1])-1
        pos_seq[pos_list.index(fixPos)]=ii.split(':')[3]
    alt_seq=''.join(pos_seq[(pos_list.index(pos-cutoff+1)):(pos_list.index(pos+cutoff+1))])
    a.append(alt_seq)
    f1.write('\t'.join(a)+'\n')
f1.close()