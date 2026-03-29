# Use MNVAnno for annotation
    grep all /data/jinww/mnv/analyse3/02data/02annotation/lnc.oneline.txt |grep -v intro > lncRNA_adjust.txt
# Perform secondary structure prediction for lncRNA 
    # Generate sequences for RNAsnp calculation 
        # Official documentation https://rth.dk/resources/rnasnp/software.php https://pubmed.ncbi.nlm.nih.gov/23630321/
    # Obtain sequences suitable for RNAsnp impact calculation
        python 001getSeq_forRNAsnp.py lncRNA_adjust.txt /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt chr_fix.fa 02RNAsnp/ & 
    # Use RNAsnp software to calculate secondary structure changes
        python 002runRNAsnp.py lncRNA_RNAsnp.txt ./ &
    # Split results into mnv_diff.res and snv_diff.res, p>0.2 means no significant change, otherwise significant
        python 003splitRes.py res_RNAsnp.txt ./
    # Calculate differences
        python 004getDiff_ms.py snv_diff.res mnv_diff.res merge_res.txt