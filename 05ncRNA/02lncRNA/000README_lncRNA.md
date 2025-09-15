# lncRNA data download
    # Download annotation information and sequence information from NONCODE
        sed -i 's/^M//g' NONCODEv6_human_hg38_lncRNA.gtf # ^M is Ctrl+V+M
    # Use getncRNAdb.py to obtain lncRNA library file
    # Already revised MNVAnnoNon.py in MNVAnno, added exon annotation
    # Soft link MNV results
        ln -s /data/jinww/mnv/analyse3/02data/01origin/hg38mnv ./
    # Get sequence data
        ln -s /data/jinww/mnv/library/human/01gtf/ensembl/hg38/chr_fix.fa ./
# Use MNVAnno for annotation
    # Only extract rows that completely fall in lnc
    grep all /data/jinww/mnv/analyse3/02data/02annotation/lnc.oneline.txt |grep -v intro > lncRNA_adjust.txt
# Perform secondary structure prediction for lncRNA [001getSeq_forRNAsnp.py for RNAsnp format, 001getSeq.py for plotting]
    # Generate sequences for RNAsnp calculation 
        # Official documentation https://rth.dk/resources/rnasnp/software.php https://pubmed.ncbi.nlm.nih.gov/23630321/
    # Obtain sequences suitable for RNAsnp impact calculation
        python 001getSeq_forRNAsnp.py lncRNA_adjust.txt /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt chr_fix.fa 02RNAsnp/ & 
    # Use RNAsnp software to calculate secondary structure changes
        python 002runRNAsnp.py lncRNA_RNAsnp.txt ./ &
    # Result interpretation
    len_seq<1000  select model1, p-value select first
        # Single snp
            SNP     W   Slen   GC   interval d_max  p-value interval r_min  p-value
            U1013C  200 3344 0.5411 975-1025 0.2432 0.0724  998-1052 0.0615 0.0932
        # Multiple snp
            SNP             W   Slen   GC   interval  d_max  p-value interval  r_min  p-value
            C9294A-U9296G   200 9605 0.4814 9261-9310 0.1951 0.0749  9268-9317 0.2345 0.1213  
    len_seq>=1000  select model2, p-value select last
        # Single snp
            SNP     w   Slen   GC   max_k d_max  p-value interval   d       p-value
            U1013C  200 3344 0.5411 994   4.3961 0.2176  994-1019   0.1265  0.1232
        # Multiple snp
            SNP             w   Slen   GC   max_k d_max  p-value interval   d       p-value
            C9294A-U9296G   200 9605 0.4814 9270  7.0487 0.0624  9270-9298  0.2463  0.0099
    # Split results into mnv_diff.res and snv_diff.res, p>0.2 means no significant change, otherwise significant
        python 003splitRes.py res_RNAsnp.txt ./
        # id,CG,p_value,flag
    # Calculate differences
        python 004getDiff_ms.py snv_diff.res mnv_diff.res merge_res.txt

# Draw pictures based on significant results
    # Extract significant MNVs first,
        python 001getSig.py ../02RNAsnp/mnv_diff.res ../01data/lncRNA_adjust.txt res_sig.mnv  # vs ref
        python 001getSig.py ../02RNAsnp/merge_res.txt ../01data/lncRNA_adjust.txt res_sig.mnv # vs SNV and ref
    # Extract sequences of significant MNVs to generate mutation before and after sequences for plotting
        python 001getSeq.py ../03RNAfold/res_sig.mnv /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt chr_fix.fa ../03RNAfold/ # This RNAFold for plotting pictures
        python 001getSeq.py ../03RNAfold_vs_SNVandRef/res_sig.mnv /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt chr_fix.fa ../03RNAfold_vs_SNVandRef/ 
        # Below to know energy
            RNAfold -p < lncRNA_MNV.seq --jobs=20 --noDP --noPS > RNAFold_lncRNA_MNV.txt
            --jobs=20: Specify the number of parallel jobs. This parameter is used to speed up prediction. Here set to 20 means run 20 jobs simultaneously.
            --noDP: Disable dynamic programming algorithm (Dynamic Programming), which can speed up prediction. Dynamic programming is used to find the optimal secondary structure, but may consume a lot of computing resources.
            --noPS: Disable output PS format scatter plot, if you don't need this format graphics, you can use this parameter to improve prediction speed.
        # Draw secondary structure diagram (fast)
            # Only draw lncRNA within 5k length, and only draw MNV and ref
            python 002lncRNA5k.py lncRNA_MNV.seq lncRNA_MNV_5k.seq
            python 002lncRNA5k.py lncRNA_ref.seq lncRNA_ref_5k.seq
            RNAfold --jobs=30 --noPS < lncRNA_MNV_5k.seq  > RNAFold_lncRNA_MNV_5k.txt &
            RNAfold --jobs=30 --noPS < lncRNA_ref_5k.seq   > RNAFold_lncRNA_ref_5k.txt &
            RNAplot -o svg < RNAFold_lncRNA_MNV.txt


