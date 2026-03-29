## Compared the MNV differences among different populations.
## Drew Venn diagrams for different populations
    python 01getMNVVenn.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv ./venn/venn_1000G.txt

## Checked the SNV Venn diagrams for different populations
    # Obtained SNV data for different populations
        python 02get1000GSNV.py /data/jinww/04reference/publicDB/1000G/hg38/ snv/
    # Counted the overlap of SNVs among different populations
        python 03getSNVVenn.py
        Obtained the probability of unique SNVs, and found the reason for the above Venn diagram results was due to the distribution of SNPs.
    # Checked the density of SNVs across different populations (all five populations identified)
        python 04getDes.py snv/AFR.snv snv/AFR_snv_density

## Extracted the MNVs from the SNVs shared across different populations and used them to draw the actual Venn diagrams for different populations
    python 05getMNVfromSameSNV.py /data/jinww/mnv/analyse3/02pop/01overlap/ /data/jinww/mnv/analyse3/02data/01origin/hg38mnv
    # Obtained common_snv and venn_1000G_common