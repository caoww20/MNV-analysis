# python getPromotor_mRNA.py en_GRCh38.108_gene.anno.txt en_GRCh38.108_gene_promotor.txt 
import sys
input_file = sys.argv[1] # /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt
output_file = sys.argv[2]

# ENST00000361624 MT      +       ENSG00000198804 MT-CO1  5904    7445    5904    7445
mydata=[]
with open(input_file) as f:
    for line in f:
        line = line.strip().split('\t')
        if line[2] == '+':
            promotor_start = int(line[5]) - 2500
            promotor_end = int(line[5])
        else:
            promotor_start = int(line[6])
            promotor_end = int(line[6]) + 2500
        
        # gene_promotor        4       47701   47943   ENST00000361624:4-47701-47943   ENSG00000198804       +
        newa=['gene_promotor',line[1],str(promotor_start),str(promotor_end),line[0]+':'+line[1]+'-'+str(promotor_start)+'-'+str(promotor_end),line[3],line[2]]
        mydata.append(newa)
# Sort mydata by the 2nd (chrom) and 3rd (start) columns
mydata.sort(key=lambda x:(x[1],int(x[2])))

with open(output_file,'w') as w:
    for newa in mydata:
        w.write('\t'.join(newa)+'\n')