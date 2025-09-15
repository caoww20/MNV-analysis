# Split aggregated MNV dataset into per-dataset files and flags
# Usage: 01splitData.py <input> <out_dir>
import sys
input=sys.argv[1]
file_url=sys.argv[2]

## Write per-dataset splits
f=open(input)
flag=open(file_url+'flag_datasets','w')
flag.write('chr\tpos\tmnvid\tref\talt\trsid\tdistance\tmnvtype\t1000G\tGTEx\tUKB20w\tUKB50w\tTCGA\tGnomAD_WGS\tGnomAD_WES\tSum\n')

GnomAD_WGS=open(file_url+'GnomAD_WGS.txt','w')
GnomAD_WGS_mnv=open(file_url+'GnomAD_WGS_mnv.txt','w')
GnomAD_WGS_multi=open(file_url+'GnomAD_WGS_mnv_multi.txt','w')

GnomAD_WES=open(file_url+'GnomAD_WES.txt','w')
GnomAD_WES_mnv=open(file_url+'GnomAD_WES_mnv.txt','w')
GnomAD_WES_multi=open(file_url+'GnomAD_WES_mnv_multi.txt','w')

G1000=open(file_url+'1000G.txt','w')
G1000_mnv=open(file_url+'1000G_mnv.txt','w')
G1000_multi=open(file_url+'1000G_mnv_multi.txt','w')

GTEx=open(file_url+'GTEx.txt','w')
GTEx_mnv=open(file_url+'GTEx_mnv.txt','w')
GTEx_multi=open(file_url+'GTEx_mnv_multi.txt','w')

UKB20w=open(file_url+'UKB20w.txt','w')
UKB20w_mnv=open(file_url+'UKB20w_mnv.txt','w')
UKB20w_multi=open(file_url+'UKB20w_mnv_multi.txt','w')

UKB50w=open(file_url+'UKB50w.txt','w')
UKB50w_mnv=open(file_url+'UKB50w_mnv.txt','w')
UKB50w_multi=open(file_url+'UKB50w_mnv_multi.txt','w')

TCGA=open(file_url+'TCGA.txt','w')
TCGA_mnv=open(file_url+'TCGA_mnv.txt','w')
TCGA_multi=open(file_url+'TCGA_mnv_multi.txt','w')
for i in f:
    a=i.strip().split('\t')
    myhead=a[:-1]
    myflag=[0,0,0,0,0,0,0]
    ismulti=len(a[2].split('.'))
    b=a[-1].split('|')
    for ii in b:
        if ii.startswith('1000G'):
            if ismulti==2:
                G1000_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                G1000.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                G1000_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')
                G1000.write('\t'.join(myhead)+'\t'+ii+'\n')
            myflag[0]=1
        elif ii.startswith('GTEx'):
            if ismulti==2:
                GTEx_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                GTEx.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                GTEx_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')
                GTEx.write('\t'.join(myhead)+'\t'+ii+'\n')
            myflag[1]=1
        elif ii.startswith('UKB20w'):
            if ismulti==2:
                UKB20w_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                UKB20w.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                UKB20w_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')
                UKB20w.write('\t'.join(myhead)+'\t'+ii+'\n')
            myflag[2]=1
        elif ii.startswith('UKB50w'):
            if ismulti==2:
                UKB50w_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                UKB50w.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                UKB50w_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')
                UKB50w.write('\t'.join(myhead)+'\t'+ii+'\n')
            myflag[3]=1
        elif ii.startswith('TCGA'):
            if ismulti==2:
                TCGA_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                TCGA.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                TCGA_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')       
                TCGA.write('\t'.join(myhead)+'\t'+ii+'\n')  
            myflag[4]=1
        elif ii.startswith('GnomAD_WGS'):
            if ismulti==2:
                GnomAD_WGS_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                GnomAD_WGS.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                GnomAD_WGS_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')
                GnomAD_WGS.write('\t'.join(myhead)+'\t'+ii+'\n')
            myflag[5]=1
        elif ii.startswith('GnomAD_WES'):
            if ismulti==2:
                GnomAD_WES_multi.write('\t'.join(myhead)+'\t'+ii+'\n')
                GnomAD_WES.write('\t'.join(myhead)+'\t'+ii+'\n')
            else:
                GnomAD_WES_mnv.write('\t'.join(myhead)+'\t'+ii+'\n')
                GnomAD_WES.write('\t'.join(myhead)+'\t'+ii+'\n')
            myflag[6]=1
    sumFlag=sum(myflag)
    myflag = [ str(x) for x in myflag ]
    flag.write('\t'.join(myhead+myflag+[str(sumFlag)])+'\n')

GnomAD_WGS.close()
GnomAD_WGS_mnv.close()
GnomAD_WGS_multi.close()
GnomAD_WES.close()
GnomAD_WES_mnv.close()
GnomAD_WES_multi.close()
G1000.close()
G1000_mnv.close()
G1000_multi.close()
GTEx.close()
GTEx_mnv.close()
GTEx_multi.close()
UKB20w.close()
UKB20w_mnv.close()
UKB20w_multi.close()
UKB50w.close()
UKB50w_mnv.close()
UKB50w_multi.close()
TCGA.close()
TCGA_mnv.close()
TCGA_multi.close()
