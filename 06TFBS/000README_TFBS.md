# new ################################################
# Some basic data
    python 001getPromotor_mRNA.py  /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt en_GRCh38.108_gene_promotor.txt 
    python 001getPromotor_mRNA.py  /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_canonical_gene.anno.txt en_GRCh38.108_canonical_gene_promotor.txt 
    python 001getPromotor_lncRNA.py /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt hg38_lncRNA_promotor.txt
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D en_GRCh38.108_canonical_gene_promotor.txt -o canical.promotor &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D en_GRCh38.108_gene_promotor.txt -o protein.promotor &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D hg38_lncRNA_promotor.txt -o lncRNA.promotor &
    python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i canical.promotor &
    python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i protein.promotor &
    python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i lncRNA.promotor &
    grep all canical.promotor.oneline.txt > canical.promotor.oneline.fix &
    grep all protein.promotor.oneline.fix.txt > protein.promotor.oneline.fix &
    grep all lncRNA.promotor.oneline.txt > lncRNA.promotor.oneline.fix & 

    python 01filter.py jaspar_hg38 jaspar_hg38_fix

# Statistics show that in protein-coding gene promoter regions.
    python 002getGeneMNVnum.py protein.promotor.oneline.fix geneMNVNum.txt
# Therefore, we collected JASPAR TF pwm matrices, used FIMO to scan TFBS before and after changes on promoter sequences, obtained xx MNV-TFBS pairs, gain has how many, loss has how many # Two methods get gain, loss numbers (take union)
    # Get union
    python ../003getTFBS.py protein.promotor.oneline.fix ../03caow/jaspar_hg38_fix jaspar jaspar_gene.txt jaspar_pair.txt &
    # Count union
    python ../004calTFBSNum.py jaspar_gene.txt jaspar_gene.num
# Average/median number of MNVs per gene affecting TFBS, xx changed MNV-TFBS pairs
    See R
# We also further evaluated the impact different from SNV, found xx% of MNV-TFBS different from SNV-TFBS pairs, and used alluvial plot to further show this change.
    python ../005match.py ../03caow/jaspar_snv jaspar_pair.txt ./jaspar_ &
    
# In addition, we also performed the same prediction for lncRNA promoter regions. 
    python ../002getGeneMNVnum.py lncRNA.promotor.oneline.fix geneMNVNum.txt
    # Get union
    python ../003getTFBS.py lncRNA.promotor.oneline.fix ../03caow/jaspar_hg38_fix jaspar jaspar_gene.txt jaspar_pair.txt &
    # Count union
    python ../004calTFBSNum.py jaspar_gene.txt jaspar_gene.num &
