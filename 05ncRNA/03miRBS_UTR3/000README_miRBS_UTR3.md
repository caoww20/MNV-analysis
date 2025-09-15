# Impact of MNVs falling on gene miRBS======================
## Get gene miRNA sequences and basic data
    # Download miRNA seed data from targetscan
        python 001getmiR.py /data/jinww/04reference/publicDB/miRNA/targetScan/hsa_miR.txt  ./ # Get all results
    # Soft link MNV results
        ln -s /data/jinww/mnv/analyse3/02data/01origin/hg38mnv ./
    # Get sequence data
        ln -s /data/jinww/mnv/library/human/01gtf/ensembl/hg38/chr_fix.fa ./
    # Get file containing gene UTR3
        ln -s /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt ./
## Extract MNVs annotated by MNVAnno falling in UTR3
    ln -s /data/jinww/mnv/analyse3/02data/02annotation/gene.oneline.txt ./
    # Only extract rows that completely fall in UTR3
    awk '$13=="UTR3_variant"' gene.oneline.txt  > UTR3_adjust.txt
## Get sequences
    # Using original full-length prediction is particularly slow, so change to ±50bp region for prediction
    python 002getSeq_50bp_UTR3.py UTR3_adjust.txt en_GRCh38.108_gene.anno.txt chr_fix.fa ./ &
## Use TargetScan for prediction
    # Prediction
        perl ../targetscan_70.pl /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_targetscan.seq tar_ref.seq res_ref.txt &
        perl ../targetscan_70.pl /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_targetscan.seq tar_SNV.seq res_SNV.txt &
        perl ../targetscan_70.pl /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_targetscan.seq tar_MNV.seq res_MNV.txt &
    # Correct position and delete non-valid miRBS
        python ../002adjust_pos.py res_ref.txt res_ref_adjust.txt &
        python ../002adjust_pos.py res_MNV.txt res_MNV_adjust.txt &
        python ../002adjust_pos.py res_SNV.txt res_SNV_adjust.txt &
    # Judge (only extract MNVs that completely fall in seed region), if ref both methods are F, but alt has one T, then T; vice versa
        python ../003getDiff_mnv.py res_ref_adjust.txt res_MNV_adjust.txt diff_MNV_vs_ref.txt &
        python ../003getDiff_snv.py res_ref_adjust.txt res_SNV_adjust.txt diff_SNV_vs_ref.txt &
    # Compared to SNV, MNV's impact on TFBS
        python ../003getDiff_ms.py  diff_SNV_vs_ref.txt diff_MNV_vs_ref.txt diff_MNV_vs_SNV.txt &
## 【Discard】Use miRmap for prediction, too slow
## Use miranda for prediction
    # Prediction
        # Strict /data/jinww/03software/miRanda-3.3a/bin/miranda miRNA_mirand.seq UTR3_mirand.seq -out miranda_res -sc 150.0 -en -30 -scale 4.0 -go -2.0 -ge -8.0 -quiet
        /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_ref.seq -out res_ref.txt -en -10 -quiet &
        /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_MNV.seq -out res_MNV.txt -en -10 -quiet &
        /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_SNV.seq -out res_SNV.txt -en -10 -quiet &
        This is a command line instruction for running a software. The meanings of each parameter are as follows:
            miranda: The software being run.
            miRNA9606.fa: The miRNA sequence file to be predicted (i.e., small RNA sequence).
            SEMA3B.fa: The target gene or target gene sequence file.
            -out All_miRNA.expressed.fa.out.0: The output file name, here named All_miRNA.expressed.fa.out.0.
            -sc 150.0: Maximum allowed score (i.e., TG value, which is the sequence heterogeneity, between 0~1000), default is 140.0, here set to 150.0.
            -en -30: Minimum free energy, the minimum free energy of sequence pairing (kcal/mol), default is 1, here set to -30. -en smaller like -50 will get fewer results. Generally, common energy thresholds range from -10 to -25
            -scale 4.0: Scale, used to adjust the relative proportion of sequence score and free energy score, default is 4.0, here set to 4.0.
            -go -2.0: Gap opening penalty parameter, default is -9.0, here set to -2.0. -go -ge smaller like -4 ge -9 will get fewer results
            -ge -8.0: Gap extending penalty parameter, default is -4.0, here set to -8.0.
            -quiet: Silent mode, does not display software running process details (warnings and errors are still displayed).
    # Extract information
        python ../002adjust_pos.py res_ref.txt res_ref_adjust.txt &
        python ../002adjust_pos.py res_SNV.txt res_SNV_adjust.txt &
        python ../002adjust_pos.py res_MNV.txt res_MNV_adjust.txt &
            # Extract key information from it
            grep '>' All_miRNA.expressed.fa.out.0 >All_miRNA.expressed.fa.out.0.miRanda.aln
            cat mresult |grep '^>>'|sed  's/^M//g'| sed 's/>>//g' >mresult.simple
            # The final file is miRNA and target gene, then input to cytoscape to draw network graph.
            miranda miRNA9606.fa SEMA3B.fa | grep '^>>'|sed  's/^M//g'| sed 's/>>//g' >miRanda_SEMA3B_final_res.txt # ^M is Ctrl+V+M

            # Result interpretation
            =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            Forward:     Score: 143.000000  Q:2 to 16  R:10 to 31 Align Len (14) (78.57%) (85.71%)

            Query:    3' uugauauGUUGGAUGAUGGAGu 5'
                                || ||| ||:|||| 
            Ref:      5' cucagguCAUCCUCCUGCCUCg 3'

            Energy:  -19.290001 kCal/Mol

            Scores for this hit:
            >hsa-let-7a-5p  rs41312668,rs1416698324|MNV01265731|ENST00000370112|+|2802|2852 143.00  -19.29  2 16    10 31   14      78.57%  85.71%

            Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions
            >>hsa-let-7a-5p rs41312668,rs1416698324|MNV01265731|ENST00000370112|+|2802|2852 143.00  -19.29  143.00  -19.29  126     22      101      10
            How to interpret the result
            This result indicates that miRanda used the input miRNA sequence and target gene sequence for prediction, and obtained a match. Specifically, the prediction results are as follows:
            A binding site for hsa-let-7a-5p miRNA was predicted, matching a region of the reference sequence, with the region starting at position 10 and ending at position 31.
            The matching situation of miRNA to the binding site is as follows: In the miRNA sequence, positions 2-16 bases form complementary pairing with positions 10-31 in the reference sequence, with some bases not completely consistent, indicating a certain error tolerance.
            The binding energy is -19.29 kcal/mol, indicating that the interaction strength between miRNA and target gene is relatively high, possibly with a certain binding relationship.
            For this prediction, miRanda output some other statistical information. Including the pairing length of 14 bases, the matching percentage of 78.57% (with miRNA sequence), and the matching percentage of 85.71% (with target gene sequence).
            In summary, this result indicates a possible miRNA-target gene interaction event, providing some quantitative information. Further analysis and verification need to be combined with more prediction and experimental data.
    # Judge (only extract MNVs that completely fall in seed region), if ref both methods are F, but alt has one T, then T; vice versa
        python ../003getDiff_mnv.py res_ref_adjust.txt res_MNV_adjust.txt diff_MNV_vs_ref.txt &
        python ../003getDiff_snv.py res_ref_adjust.txt res_SNV_adjust.txt diff_SNV_vs_ref.txt &
    # Compared to SNV, MNV's impact on TFBS
        python ../003getDiff_ms.py  diff_SNV_vs_ref.txt diff_MNV_vs_ref.txt diff_MNV_vs_SNV.txt &

# Impact of MNVs falling on lncRNA======================
## Extract MNVs annotated by MNVAnno falling in lncRNA exon
    ln -s  /data/jinww/mnv/analyse3/02data/02annotation/lnc.oneline.txt  lncRNA_adjust.txt
## Get lncRNA sequences
    python 002getSeq_50bp_lncRNA.py lncRNA_adjust.txt /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt chr_fix.fa ./ &
## Use TargetScan for prediction
    # Prediction
        perl ../targetscan_70.pl /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_targetscan.seq tar_ref.seq res_ref.txt &
        perl ../targetscan_70.pl /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_targetscan.seq tar_SNV.seq res_SNV.txt &
        perl ../targetscan_70.pl /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_targetscan.seq tar_MNV.seq res_MNV.txt &
## Use miranda for prediction
    /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_ref.seq -out res_ref.txt -en -10 -quiet &
    /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_MNV.seq -out res_MNV.txt -en -10 -quiet &
    /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_SNV.seq -out res_SNV.txt -en -10 -quiet &

# Merge results from two methods
## First level merge: two methods consistent as consistent, if both loss and gain then mark as unclear, only one then record one
    python 001mergeRes.py ../02miranda/UTR3/diff_MNV_vs_ref.txt ../03targetscan/UTR3/diff_MNV_vs_ref.txt UTR3/diff_MNV_vs_ref.txt &
    python 001mergeRes.py ../02miranda/UTR3/diff_SNV_vs_ref.txt ../03targetscan/UTR3/diff_SNV_vs_ref.txt UTR3/diff_SNV_vs_ref.txt &
    python 001mergeRes.py ../02miranda/lncRNA/diff_MNV_vs_ref.txt ../03targetscan/lncRNA/diff_MNV_vs_ref.txt lncRNA/diff_MNV_vs_ref.txt &
    python 001mergeRes.py ../02miranda/lncRNA/diff_SNV_vs_ref.txt ../03targetscan/lncRNA/diff_SNV_vs_ref.txt lncRNA/diff_SNV_vs_ref.txt &
## Compared to SNV, MNV's impact on TFBS, containing miss results not considered, for SNV if one none one gain then gain, different then unclear (also not considered)
    python 002getDiff.py  UTR3/diff_SNV_vs_ref.txt UTR3/diff_MNV_vs_ref.txt UTR3/diff_MNV_vs_SNV.txt &
    python 002getDiff.py  lncRNA/diff_SNV_vs_ref.txt lncRNA/diff_MNV_vs_ref.txt lncRNA/diff_MNV_vs_SNV.txt &
    # grep gain diff_MNV_vs_SNV.txt |grep loss  theoretically no results, if have then problem!!

# Get some other results
    ## Get number of MNVs on UTR3
    python 003getGeneNum_UTR3.py /data/jinww/mnv/analyse3/02data/02annotation/gene.oneline.txt res/UTR3_MNVNum.txt 
    ## Get number of different MNVs on UTR3 (different from ref)
    python 003getGeneNum_UTR3_vs_ref.py /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt UTR3/diff_MNV_vs_ref.txt res/UTR3_MNVNum_vs_ref.txt
    ## See number of different MNVs on UTR3 (MNV different from SNV)
    python 003getGeneNum_UTR3_vs_SNV.py /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt UTR3/diff_MNV_vs_SNV.txt res/UTR3_MNVNum_vs_SNV.txt
    ## Draw alluvial
    python 004alluvial.py UTR3/diff_MNV_vs_SNV.txt res/UTR3_alluvial.txt

    ## Get number of MNVs on lncRNA (only all, and only exon)
    python 005getGeneNum_lncRNA.py /data/jinww/mnv/analyse3/02data/02annotation/lnc.oneline.txt res/lncRNA_MNVNum.txt 
    ## Get number of different MNVs on UTR3 (different from ref)
    python 005getGeneNum_lncRNA_vs_ref.py /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt lncRNA/diff_MNV_vs_ref.txt res/lncRNA_MNVNum_vs_ref.txt
    ## See number of different MNVs on UTR3 (MNV different from SNV)
    python 005getGeneNum_lncRNA_vs_SNV.py /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt lncRNA/diff_MNV_vs_SNV.txt res/lncRNA_MNVNum_vs_SNV.txt
    ## Draw alluvial
    python 004alluvial.py lncRNA/diff_MNV_vs_SNV.txt res/lncRNA_alluvial.txt



# Impact of MNVs falling on mature miRNA on miRBs======================
## Get falling on mature miRNA
    ln -s /data/jinww/mnv/analyse3/05ncRNA/01pre_miRNA/01data/mature_miRNA_adjust.txt ../05altMIR/
## Generate changed miRNA sequences, only care about seed region
    python 002getSeq_miRNA.py mature_miRNA_adjust.txt chr_fix.fa ./
## Get complete UTR3 sequences of genes and lncRNA sequences (here only need to extract ref)
    # Extract full-length UTR3 then transfer
        python 002getSeq_total_UTR3.py ../01data/en_GRCh38.108_gene.anno.txt ../01data/chr_fix.fa ./ &
    # Extract full-length lncRNA, must run above first because name same, then transfer location
        python 002getSeq_total_lncRNA.py /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt ../01data/chr_fix.fa ./ &
## Use TargetScan for prediction
    # tar UTR3/lncRNA
        perl ../targetscan_70.pl ../../miR_tar_ref.seq tar_ref.seq res_ref.txt &
        perl ../targetscan_70.pl ../../miR_tar_MNV.seq tar_ref.seq res_MNV.txt &
    # Get difference
        python ../003getDiff.py res_ref.txt res_MNV.txt diff_MNV_vs_ref.txt &
## Use miRmap for prediction
    # miranda UTR3/lncRNA
        /data/jinww/03software/miRanda-3.3a/bin/miranda ../../miR_miranda_ref.seq miranda_ref.seq -out res_ref.txt -en -10 -quiet &
        /data/jinww/03software/miRanda-3.3a/bin/miranda ../../miR_miranda_MNV.seq miranda_ref.seq -out res_MNV.txt -en -10 -quiet &
    # Extract data
        python ../002adjust_pos.py res_ref.txt res_ref_adjust.txt &
        python ../002adjust_pos.py res_MNV.txt res_MNV_adjust.txt &
    # Get difference
        python ../003getDiff.py res_ref_adjust.txt res_MNV_adjust.txt diff_MNV_vs_ref.txt &
## Merge results, only retain consistent results
    python 003mergeRes.py miranda/UTR3/diff_MNV_vs_ref.txt targetscan/UTR3/diff_MNV_vs_ref.txt merge/UTR3/diff_MNV_vs_ref.txt
    python 003mergeRes.py miranda/lncRNA/diff_MNV_vs_ref.txt targetscan/lncRNA/diff_MNV_vs_ref.txt merge/lncRNA/diff_MNV_vs_ref.txt

# Draw pictures
    # Secondary structure
        premiRNA difference score barplot demo
        lncRNA difference grade barplot demo
    # MNV-miRBS impact
        Number of MNVs on UTR3
        Number of gains and losses, proportion different from SNV, pie chart
        Impact map
    # Seed region MNV impact
        17 MNVs after change can bind or lose how many gene bindings

# Some tests
## Test if extraction has problem
    ENST00000308647 '+'
    ENST00000217185 '-'
    awk '$3=="-"' /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt |awk '{print $1,$3,$17}'|less -S
    No problem
## Test if adding 50 has problem
    rs766256036,rs757730549|MNV06425718_MI0016596_-
    awk '$3=="-"' /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt |awk '{print $1,$3,$17}'|less -S
    No problem





