# python 001getPopSample.py /data/jinww/04reference/publicDB/1000G/igsr-1000_genomes_30x_on_grch38.tsv  /data/jinww/mnv/analyse2/05population/02group/01origin/sample/
import sys

file=sys.argv[1] #  file_common='/data/jinww/04reference/publicDB/1000G/igsr-1000_genomes_30x_on_grch38.tsv'
res=sys.argv[2] # file_origin='/data/jinww/mnv/analyse2/05population/02group/01origin/sample/'

pop1={}
pop2={}

with open(file) as f:
    next(f)
    for i in f:
        a=i.strip().split('\t')
        if a[3] not in pop1:
            pop1[a[3]]=[a[0]]
        else:
            pop1[a[3]].append(a[0])
        if a[5] not in pop2:
            pop2[a[5]]=[a[0]]
        else:
            pop2[a[5]].append(a[0])

for i in pop1:
    with open(res+i+'.sample','w') as w:
        for j in pop1[i]:
            w.write(j+'\n')
    with open(res+'pop1.sample','a+') as w:
        for j in pop1[i]:
            if i!='IBS,MSL':
                w.write(i+'\t'+j+'\n')    
        

for i in pop2:
    with open(res+i+'.sample','w') as w:
        for j in pop2[i]:
            w.write(j+'\n')
    with open(res+'pop2.sample','a+') as w:
        for j in pop2[i]:
            if i!='EUR,AFR':
                w.write(i+'\t'+j+'\n')

