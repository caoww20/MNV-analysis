import os,sys
url=sys.argv[1] # '/data/jinww/04reference/publicDB/1000G/hg38/'
res_url=sys.argv[2] # './'

file_list=os.listdir(url)
file_list = [s for s in file_list if 'vcf' in s]

AFR=open(res_url+'AFR.snv','w')
AMR=open(res_url+'AMR.snv','w')
EUR=open(res_url+'EUR.snv','w')
EAS=open(res_url+'EAS.snv','w')
SAS=open(res_url+'SAS.snv','w')
for vcf in file_list:
    f=open(url+vcf)
    for i in f:
        if i[0]!='#':
            a=i[0:2000].split('\t')
            b=a[7].split(';')
            if len(a[3])==1 and len(a[4])==1:
                for ii in b:
                    if ii.startswith('AC_AFR=') and ii.split('=')[1]!='0':
                        AFR.write('\t'.join(['AFR',a[2],a[0].replace('chr',''),a[1]])+'\n')
                    if ii.startswith('AC_AMR=') and ii.split('=')[1]!='0':
                        AMR.write('\t'.join(['AMR',a[2],a[0].replace('chr',''),a[1]])+'\n')
                    if ii.startswith('AC_EUR=') and ii.split('=')[1]!='0':
                        EUR.write('\t'.join(['EUR',a[2],a[0].replace('chr',''),a[1]])+'\n')
                    if ii.startswith('AC_EAS=') and ii.split('=')[1]!='0':
                        EAS.write('\t'.join(['EAS',a[2],a[0].replace('chr',''),a[1]])+'\n')
                    if ii.startswith('AC_SAS=') and ii.split('=')[1]!='0':
                        SAS.write('\t'.join(['SAS',a[2],a[0].replace('chr',''),a[1]])+'\n')
    f.close()
AFR.close()
AMR.close()
EUR.close()
EAS.close()
SAS.close()