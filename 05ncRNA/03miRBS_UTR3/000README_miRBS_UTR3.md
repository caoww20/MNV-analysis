# Impact of MNVs falling on gene miRBS======================
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
## Use miranda for prediction
    # Prediction
        /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_ref.seq -out res_ref.txt -en -10 -quiet &
        /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_MNV.seq -out res_MNV.txt -en -10 -quiet &
        /data/jinww/03software/miRanda-3.3a/bin/miranda /data/jinww/mnv/analyse3/05ncRNA/03miRBS_UTR3/01data/miR_mirand.seq miranda_SNV.seq -out res_SNV.txt -en -10 -quiet &
    # Extract information
        python ../002adjust_pos.py res_ref.txt res_ref_adjust.txt &
        python ../002adjust_pos.py res_SNV.txt res_SNV_adjust.txt &
        python ../002adjust_pos.py res_MNV.txt res_MNV_adjust.txt &
    # Judge
        python ../003getDiff_mnv.py res_ref_adjust.txt res_MNV_adjust.txt diff_MNV_vs_ref.txt &
        python ../003getDiff_snv.py res_ref_adjust.txt res_SNV_adjust.txt diff_SNV_vs_ref.txt &
    # Compared to SNV, MNV's impact
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
## Compared to SNV, MNV's impact
    python 002getDiff.py  UTR3/diff_SNV_vs_ref.txt UTR3/diff_MNV_vs_ref.txt UTR3/diff_MNV_vs_SNV.txt &
    python 002getDiff.py  lncRNA/diff_SNV_vs_ref.txt lncRNA/diff_MNV_vs_ref.txt lncRNA/diff_MNV_vs_SNV.txt &







