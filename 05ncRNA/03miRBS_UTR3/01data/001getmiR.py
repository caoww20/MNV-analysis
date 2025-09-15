import sys

miR = sys.argv[1] # hsa_miR.txt
output_file = sys.argv[2] # utr3/

####################################################################
# let-7-3p/miR-3596d/98-3p        UAUACAA 10116   rno-let-7a-1-3p CUAUACAAUCUACUGUCUUUCC  -1      MIMAT0017085

# targetscan
# hsa-miR-9500    AGGGAAG 9606

# miranda
# >hsa-miR-9500
# AAGGGAAGAUGGUGACCAC

# mirmap
# hsa-let-7a-5p   UGAGGUAGUAGGUUGUAUAGUU  9606

# Obtain 5'utr and cds sequences
with open(miR) as f, open(output_file+'miR_targetscan.seq','w') as result1, open(output_file+'miR_mirand.seq','w') as result2:
    for i in f:
        if 'Accession' not in i:
            line=i.strip('\n').split('\t')
            result1.write('\t'.join([line[3],line[1],'9606'])+'\n')
            result2.write('>'+line[3]+'\n')
            result2.write(line[4]+'\n')
####################################################################