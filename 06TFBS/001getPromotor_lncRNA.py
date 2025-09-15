# python getPromotor_mRNA.py hg38_lncRNA.txt en_GRCh38.108_gene_promotor.txt 
import sys
input_file = sys.argv[1] # /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt
output_file = sys.argv[2]

# ENST00000361624 MT      +       ENSG00000198804 MT-CO1  5904    7445    5904    7445
# lncRNA  1       11869   14409   +       NONHSAT148173.1 .       NONHSAG000001.2 1-11869-12227,2-12613-12721,3-13221-14409       NONCODE
mydata=[]
with open(input_file) as f:
    for line in f:
        line = line.strip().split('\t')
        if line[4] == '+':
            promotor_start = int(line[2]) - 2500
            promotor_end = int(line[2])
        else:
            promotor_start = int(line[3])
            promotor_end = int(line[3]) + 2500
        
        # gene_promotor        4       47701   47943   ENST00000361624:4-47701-47943   ENSG00000198804       +
        newa=['lncRNA_promotor',line[1],str(promotor_start),str(promotor_end),line[5]+':'+line[1]+'-'+str(promotor_start)+'-'+str(promotor_end),line[7],line[4]]
        mydata.append(newa)
# Sort mydata by the 2nd (chrom) and 3rd (start) columns
mydata.sort(key=lambda x:(x[1],int(x[2])))

with open(output_file,'w') as w:
    for newa in mydata:
        w.write('\t'.join(newa)+'\n')