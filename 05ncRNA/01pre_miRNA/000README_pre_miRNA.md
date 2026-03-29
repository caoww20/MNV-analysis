# Use MNVAnno for annotation
    # pre miRNA
        grep all /data/jinww/mnv/analyse3/02data/02annotation/pre_miRNA.oneline.txt > pre_miRNA_adjust.txt # Only extract rows that completely fall in pre miRNA
    # mature miRNA
        grep all /data/jinww/mnv/analyse3/02data/02annotation/mature_miRNA.oneline.txt > mature_miRNA_adjust.txt
    
# Perform secondary structure prediction for pre-miRNA
    # Generate sequences before and after mutation 
        python 001getSeq.py 01data/pre_miRNA_adjust.txt 01data/chr_fix.fa 02RNAFold/
    # Use RNAFold-ViennaRNA software to calculate minimum free energy of sequences before and after mutation
        https://www.tbi.univie.ac.at/RNA/ 
        https://www.tbi.univie.ac.at/RNA/documentation.html
        https://www.tbi.univie.ac.at/RNA/RNAplot.1.html
        # Install software            
            conda install -c bioconda viennarna (MNVAnno2)
        # Run command
            RNAfold -p < preMIR_ref.seq > ref.res & # This will generate dp.ps and ss.ps to get a MFE and a ps format file
            RNAfold -p < preMIR_SNV.seq > SNV.res &
            RNAfold -p < preMIR_MNV.seq > MNV.res &

# Extract important information from RNAfold results
    python 002getMFE.py 02RNAFold/mnv/mnv.res 02RNAFold/mnv/mnv_fix.res
    python 002getMFE.py 02RNAFold/snv/snv.res 02RNAFold/snv/snv_fix.res
    python 002getMFE.py 02RNAFold/ref/ref.res 02RNAFold/ref/ref_fix.res

# Calculate differences
    python 003getDiff_mnv.py 02RNAFold/ref/ref_fix.res 02RNAFold/mnv/mnv_fix.res 02RNAFold/mnv/mnv_diff.res 
    python 003getDiff_snv.py 02RNAFold/ref/ref_fix.res 02RNAFold/snv/snv_fix.res 02RNAFold/snv/snv_diff.res
    python 003getDiff_ms.py 02RNAFold/snv/snv_diff.res 02RNAFold/mnv/mnv_diff.res 02RNAFold/merge_res.txt