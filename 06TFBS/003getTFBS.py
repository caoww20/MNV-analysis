# python 003getTFBS.py protein.promotor.oneline.fix ../03caow/jaspar_hg38_fix jaspar jaspar_gene.txt jaspar_pair.txt
import sys
input1 = sys.argv[1] # input1='protein.promotor.oneline.fix'
input2 = sys.argv[2] # input2='../03caow/hocomoco_hg38_fix'
flag = sys.argv[3] # flag='jaspar'
output_file = sys.argv[4] # output_file='jaspar_gene.txt'
output_file2 = sys.argv[5] # output_file='jaspar_pair.txt'

chr_list={}
for i in range(1,23):
    chr_list[str(i)]={}
chr_list['X']={}


with open(input1) as f:
    for line in f:
        line = line.strip().split('\t')
        info=line[-1].split(' ')
        geneid=info[6]
        mnvid=line[2]
        chrid=line[0]
        if geneid in chr_list[chrid]:
            chr_list[chrid][geneid][0].append(mnvid)
        else:
            # Structure: [mnv_ids, diff_mnv_ids, diff_tfbs_ids, gain_ids, loss_ids, strand]
            chr_list[chrid][geneid]=[[mnvid],[],[],[],[],info[-1]]

for i in chr_list:
    for j in chr_list[i]:
        chr_list[i][j][0]=list(set(chr_list[i][j][0]))

# Store lines where TFBS overlap with MNVs
mydata=[]

with open(input2) as f:
    next(f)
    for line in f:
        i=line.strip().split('\t')
        chrid=i[2].replace('chr','')
        if flag=='jaspar':
            tf=i[1]+'_'
        else:
            tf=''
        mnvid=i[4]
        strand=i[8]
        mytype=i[-2]
        motif=i[0]
        TF=i[1]
        for n in chr_list[chrid]:
            # Check if on the same strand
            if strand==chr_list[chrid][n][-1]:
                # Check if the MNV exists for this gene
                if mnvid in chr_list[chrid][n][0]:
                    chr_list[chrid][n][1].append(mnvid)
                    chr_list[chrid][n][2].append(tf+motif+'|'+mnvid)
                    mydata.append(line)
                    if mytype=='Gain':
                        chr_list[chrid][n][3].append(tf+motif+'|'+mnvid)
                    else:
                        chr_list[chrid][n][4].append(tf+motif+'|'+mnvid)

# Output table: gene -> aggregated info
with open(output_file,'w') as w:
    for chrid in chr_list:
        for i in chr_list[chrid]:
            # Counts: total MNVs, diff MNVs, diff TFBSs, gains, losses, strand
            # gene structure: [[mnvid], [], [], [], [], strand]
            newa=[i,','.join(list(set(chr_list[chrid][i][0]))),','.join(list(set(chr_list[chrid][i][1]))),','.join(list(set(chr_list[chrid][i][2]))),
            ','.join(list(set(chr_list[chrid][i][3]))),','.join(list(set(chr_list[chrid][i][4]))),chr_list[chrid][i][5]]
            w.write('\t'.join(newa)+'\n')
# Output raw matched TFBS lines
mydata=set(mydata)
with open(output_file2,'w') as w:
    for i in mydata:
        w.write(i)
