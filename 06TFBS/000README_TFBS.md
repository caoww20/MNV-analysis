# Calculate overlap of 2 TF libraries
    See 000plot.R
# Calculate similarity of motifs in two libraries, how many with similarity 90% or above, p-value 0.01
# Average number of unique MNVs per gene promoter region (unique MNVs / deduplicated upstream length), calculated separately for all transcripts, protein-coding, and lncRNA (select MNV AC >1)
# Distribution of number of MNVs per gene promoter region carrying MNVs, e.g., how many genes carry 1 MNV, 2 MNVs... (select MNV AC >1, here refers to genes)
# FIMO identification results statistics, two methods two changes (gain, loss bar chart)
# How many MNVs (xx%) cause TFBS changes
# How many MNV-TFBS pairs (xx%) cause changes different from SNV-TFBS
# Draw alluvial plot: To understand differences between SNV and MNV, we first get a paired data of SNV-TFBS and MNV-TFBS 3*3 number changes, then draw the alluvial plot of SNV-TFBS -> MNV-TFBS
    python 001match.py snv.txt mnv.txt ./
# ±2.5kb divided into blocks, 100bp per unit
    # CPG curve
    # Curve of number of different MNV-TFBS
    # Boxplot of p-value curve
        First draw a boxplot curve of p-values for one TF
    # Advanced version:
        Get the median p-value for each TF in this region
        Then draw boxplot curves for all TFs
# demo
# MNVQTL enriched in cancer



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

    # In caowen's data, there may be redundant 5+ joints, discard them
    python 01filter.py hocomoco_hg38 hocomoco_hg38_fix
    python 01filter.py jaspar_hg38 jaspar_hg38_fix

# It is well known that TF binding in promoter regions can affect gene transcription, and MNVs can indeed change the sequence of binding regions to affect TF binding.
## Statistics show that in protein-coding gene promoter regions, xx (xx%) contain at least 1 MNV, xx% contain at least 5 or more MNV changes. (I guess there are no genes without promoter MNVs)
# Statistics show that in protein-coding gene promoter regions, the average number of MNVs is xx, xx% contain at least 5 or more MNV changes.
    python 002getGeneMNVnum.py protein.promotor.oneline.fix geneMNVNum.txt
# Therefore, we collected JASPAR and HOCOMOCO TF pwm matrices, used FIMO to scan TFBS before and after changes on promoter sequences, obtained xx MNV-TFBS pairs, gain has how many, loss has how many # Two methods get gain, loss numbers (take union)
    # Get union
    python ../003getTFBS.py protein.promotor.oneline.fix ../03caow/hocomoco_hg38_fix hocomoco hocomoco_gene.txt hocomoco_pair.txt &
    python ../003getTFBS.py protein.promotor.oneline.fix ../03caow/jaspar_hg38_fix jaspar jaspar_gene.txt jaspar_pair.txt &
    # Count union
    python ../004calTFBSNum.py hocomoco_gene.txt hocomoco_gene.num
    python ../004calTFBSNum.py jaspar_gene.txt jaspar_gene.num
    # Get gain loss numbers for each method, and merged gain loss numbers, total MNV-TFBS numbers (whether to take union or add directly?)
    # No deduplication, because motif differences between two methods are large
    type    gain    loss    total
    hoco    345607  411634   757241
    jaspar  369173  421359   790532
# Average/median number of MNVs per gene affecting TFBS, xx changed MNV-TFBS pairs
    See R
# We also further evaluated the impact different from SNV, found xx% of MNV-TFBS different from SNV-TFBS pairs, and used alluvial plot to further show this change.
    python ../005match.py ../03caow/jaspar_snv jaspar_pair.txt ./jaspar_ &
    python ../005match.py ../03caow/hocomoco_snv hocomoco_pair.txt ./hocomoco_ &
# We compared TFBS distribution at different distances, found that the closer to TSS, the more MNV numbers, which is consistent with CPG distribution, reinforcing the first mechanism of MNV generation.
# Similarly, we also assessed the impact of MNV on TFBS at different distances, found that the closer to TSS, a single MNV can cause more TFBS differences, which highlights the importance of MNVs close to TSS

# In addition, we also performed the same prediction for lncRNA promoter regions. Obtained xx MNV-TFBS pairs, gain has how many, loss has how many.
# lncRNA promoter regions MNV numbers average xx, with protein-coding MNV numbers not significantly different, but MNV affecting TFBS numbers slightly lower than protein, MNV impact on protein-coding may be higher than lncRNA 【Boxplot】
    python ../002getGeneMNVnum.py lncRNA.promotor.oneline.fix geneMNVNum.txt
    # Get union
    python ../003getTFBS.py lncRNA.promotor.oneline.fix ../03caow/hocomoco_hg38_fix hocomoco hocomoco_gene.txt hocomoco_pair.txt &
    python ../003getTFBS.py lncRNA.promotor.oneline.fix ../03caow/jaspar_hg38_fix jaspar jaspar_gene.txt jaspar_pair.txt &
    # Count union
    python ../004calTFBSNum.py hocomoco_gene.txt hocomoco_gene.num &
    python ../004calTFBSNum.py jaspar_gene.txt jaspar_gene.num &
    # Get gain loss numbers for each method, and merged gain loss numbers, total MNV-TFBS numbers (whether to take union or add directly?)
    # No deduplication, because motif differences between two methods are large
    type    gain    loss    total
    hoco    1032408 1218150   2250558     
    jaspar  1150383 1290093   2440476 


# In summary, we performed systematic identification of TFBS in protein-coding and lncRNA promoter regions, and confirmed through visualization that these base changes affect TFBS binding.
# A total of xx identified, comprehensive annotation of promoter region MNVs will help research on what.