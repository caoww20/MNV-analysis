# python 014getcardiomyopathyMNV.py dilated_cardiomyopathy.gwas /data/jinww/mnv/analyse3/02data/02datasets/1000G.txt overlap.MNV
import sys
inputfile1 = sys.argv[1] # inputfile1='/data/jinww/mnv/analyse3/02data/02datasets/1000G.txt'
inputfile2 = sys.argv[2] # inputfile2='dilated_cardiomyopathy.gwas'
outputfile = sys.argv[3] # outputfile='overlap_dilated.MNV'

snv=[]
with open(inputfile2) as f:
    for line in f:
        a=line.strip().split('\t')[21]
        snv.append(a)
snv=list(set(snv))

snv2=[]
with open(inputfile1) as f, open(outputfile,'w') as out:
    for line in f:
        a=line.strip().split('\t')
        rsid=a[5].split(',')
        for ii in rsid:
            if ii in snv:
                snv2.append(ii)
                out.write(line)
                break

snv2=list(set(snv2))
# print(len(snv2))
with open(inputfile2) as f, open(inputfile2+'.filter','w') as out:
    for line in f:
        for ii in snv2:
            if ii in line:
                out.write(line)
                break