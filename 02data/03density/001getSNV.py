## 获得SNV信息
import os,sys,gzip
file=sys.argv[1]
flag=sys.argv[2]
res=sys.argv[3]
def openfile(filename, mode="r"):
    if filename[-3:] == '.gz':
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
if flag=='1000G':
    file_list=os.listdir(file)
    file_list = [s for s in file_list if 'chrX' not in s]
    res=open(res,'w')
    for vcf in file_list:
        f=open(file+vcf)
        for i in f:
            if i[0]!='#':
                a=i[0:2000].split('\t')
                chr=a[0]
                pos=a[1]
                snv=a[2]
                if len(a[3])==1 and len(a[4])==1:
                    res.write('\t'.join([snv,chr.replace('chr',''),pos])+'\n')
        f.close()
    res.close()
elif flag=='GTEx':
    file_list=os.listdir(file)
    res=open(res,'w')
    for vcf in file_list:
        f=open(file+vcf)
        for i in f:
            if i[0]!='#':
                a=i[0:2000].split('\t')
                chr=a[0]
                pos=a[1]
                snv=a[2]
                if len(a[3])==1 and len(a[4])==1:
                    res.write('\t'.join([snv,chr.replace('chr',''),pos]).replace('chr','').replace('_b38','').replace('_',':')+'\n')
        f.close()
    res.close()    
elif flag=='UKB20w':
    file_list=os.listdir(file)
    file_list = [s for s in file_list if '.vcf.gz' in s and '_X_' not in s and 'reject' not in s]
    res=open(res,'w')
    for vcf in file_list:
        with openfile(file+vcf) as f:
            for i in f:
                i = i[0:3000].decode('utf-8')
                if i[0]!='#':
                    a=i[0:2000].split('\t')
                    chr=a[0].replace('chr','')
                    pos=a[1]
                    snv=a[2]
                    if len(a[3])==1 and len(a[4])==1:
                        res.write('\t'.join([snv,chr,pos])+'\n')
    res.close()        
elif flag=='UKB50w':
    file_list=os.listdir(file)
    file_list = [s+'/ukb_imp_'+s+'.vcf.gz' for s in file_list if 'chr' in s and 'XY' not in s and 'chr14' not in s]
    file_list.append('chr14/ukb_imp_chr14_filter.vcf.gz')
    res=open(res,'w')
    for vcf in file_list:
        with openfile(file+vcf) as f:
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
elif flag=='TCGA':
    file_list=os.listdir(file)
    res=open(res,'w')
    for vcf in file_list:
        f=open(file+vcf)
        for i in f:
            if i[0]!='#':
                a=i[0:2000].split('\t')
                chr=a[0]
                pos=a[1]
                snv=a[2]
                if len(a[3])==1 and len(a[4])==1:
                    res.write('\t'.join([snv,chr,pos])+'\n')
        f.close()
    res.close()    

