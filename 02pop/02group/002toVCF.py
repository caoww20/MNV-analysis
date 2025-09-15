import sys

VCF_head=sys.argv[1] #  file_common='01origin/VCF.head'
file=sys.argv[2] # file_origin='02adjust/commonMNV.genotype'
res=sys.argv[3] # res='02adjust/commonMNV.vcf'

with open(VCF_head,'r') as f, open(res,'w') as f2:
    for line in f:
        f2.write(line)

with open(file,'r') as f, open(res,'a+') as f2:
    first_line = f.readline()
    line=first_line.strip().split('\t')
    newline=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+line[7:]
    f2.write('\t'.join(newline)+'\n')

for chrid in list(range(1, 23)):
    pos_list=[]
    with open(file,'r') as f, open(res,'a+') as f2:
        next(f)
        for line in f:
            mnv_chr=int(line[0:10].split('\t')[0])
            if chrid==mnv_chr:
                line=line.strip().split('\t')
                pos=line[1].split(',')[0]
                if pos not in pos_list:
                    pos_list.append(pos)
                    newline=[line[0],pos,line[2],line[3].split(',')[0],line[4].split(',')[0],'.','PASS',line[6],'GT']+line[7:]
                else:
                    for j in range(1,10001): # 这里必须给大一些，不然有些点就输出不出来
                        new_pos=str(int(pos)+j)
                        if new_pos not in pos_list:
                            pos_list.append(new_pos)
                            newline=[line[0],new_pos,line[2],line[3].split(',')[0],line[4].split(',')[0],'.','PASS',line[6],'GT']+line[7:]
                            break
                # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
                # 1 85  mnv1    A   T   .   PASS    .   GT
                f2.write('\t'.join(newline)+'\n')