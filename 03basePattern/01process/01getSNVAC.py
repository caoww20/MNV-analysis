
import sys
filename=sys.argv[1] # /data/jinww/04reference/publicDB/
resname=sys.argv[2] # /data/jinww/mnv/analyse2/04basePattern/01process

## Extract SNV AC/AF from VCF (IGSR)
f1=open(resname+'1000G_snv.txt','w')
for chr in list(range(1,23))+['X']:
    ATCG=['A','T','C','G']
    with open(filename+'1000G/hg38/IGSR.chr'+str(chr)+'.vcf') as f:
        for i in f:
            # chr22	10519265	22:10519265:CA:C	CA	C	.	.	AC=2;AF=0.000312305;
            if i[0]!='#':
                a=i[0:1000].split('\t')
                if a[3] in ATCG and a[4] in ATCG:
                    id=':'.join([a[0].replace('chr',''),a[1],a[3],a[4]])
                    ac=a[7].split(';')[0].replace('AC=','')
                    af=a[7].split(';')[1].replace('AF=','')
                    f1.write('\t'.join([id,ac,af])+'\n')
f1.close()


## Extract SNV AC/AF from VCF (GTEx)
f1=open(resname+'GTEx_snv.txt','w')
for chr in list(range(1,23))+['X']:
    ATCG=['A','T','C','G']
    with open(filename+'GTEx/hg38/GTEx.chr'+str(chr)+'.vcf') as f:    
        for i in f:
            # 1       13526   chr1_13526_C_T_b38      C       T       .       PASS    AN=1676;AF=0.000596659;AC=1
            if i[0]!='#':
                a=i[0:1000].split('\t')
                if a[3] in ATCG and a[4] in ATCG:
                    id=':'.join([a[0],a[1],a[3],a[4]])
                    ac=a[7].split(';')[2].replace('AC=','')
                    af=a[7].split(';')[1].replace('AF=','')
                    f1.write('\t'.join([id,ac,af])+'\n')
f1.close()