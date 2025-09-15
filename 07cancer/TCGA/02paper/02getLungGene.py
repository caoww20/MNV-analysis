# python 02getLungGene.py gene.txt abstract.txt lung_gene.txt
import sys
input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]

gene={}
with open(input_file1,"r") as file:
    for line in file:
        gene[line.strip()]=0

with open(input_file2,"r") as file:
    for line in file:
        line=line.upper()
        for key in gene:
            if key in line:
                gene[key]=gene[key]+1

with open(output_file,"w") as file:
    for key in gene:
        file.write(key+'\t'+str(gene[key])+"\n")