# miRNA data download
    # Soft link MNV results
        ln -s /data/jinww/mnv/analyse3/02data/01origin/hg38mnv ./
    # Get sequence data
        ln -s /data/jinww/mnv/library/human/01gtf/ensembl/hg38/chr_fix.fa ./
# Use MNVAnno for annotation
    # pre miRNA
        grep all /data/jinww/mnv/analyse3/02data/02annotation/pre_miRNA.oneline.txt > pre_miRNA_adjust.txt # Only extract rows that completely fall in pre miRNA
    # mature miRNA
        grep all /data/jinww/mnv/analyse3/02data/02annotation/mature_miRNA.oneline.txt > mature_miRNA_adjust.txt
    
# Perform secondary structure prediction for pre-miRNA
    # Generate sequences before and after mutation # Here use _ instead of |, because considering the need for generating pictures later, | is not suitable as filename, like touch a|a.txt is not good, but a_a.txt is fine, the reason not to use - is because there is strand
        python 001getSeq.py 01data/pre_miRNA_adjust.txt 01data/chr_fix.fa 02RNAFold/
    # Use RNAFold-ViennaRNA software to calculate minimum free energy of sequences before and after mutation
        https://www.tbi.univie.ac.at/RNA/ 
        https://www.tbi.univie.ac.at/RNA/documentation.html
        https://www.tbi.univie.ac.at/RNA/RNAplot.1.html
        # Install software            
            conda install -c bioconda viennarna (MNVAnno2)
        # Run command
            # RNAfold < test.fa > test.res # This can only generate one ps, the drawing is ugly, so use the following
            RNAfold -p < preMIR_ref.seq > ref.res & # This will generate dp.ps and ss.ps to get a MFE and a ps format file
            RNAfold -p < preMIR_SNV.seq > SNV.res &
            RNAfold -p < preMIR_MNV.seq > MNV.res &
    # Draw secondary structure diagram
    # Use RNAplot to draw pictures, but ugly
        RNAplot -o svg < test.res
    # Use magick to draw pictures
        /data/jinww/03software/ViennaRNA-2.5.0/src/Utils/relplot.pl MNV01003213_ss.ps MNV01003213_dp.ps > MNV01003213_rss.ps # Reliability information annotated with color on RNA secondary structure diagram (beauty)
        # Convert postscript file to pdf (also can use ghostscript to convert)
        magick convert -density 300 5S_rss.ps 5S_rss.pdf # Generate secondary structure pdf diagram
        magick convert -density 300 5S_rss.pdf 5S_rss.png # Generate secondary structure png diagram

# Extract important information from RNAfold results
    >rs778319134,rs757730549_MNV06425718_MI0016596_-
    CUGGCUUCCAAAGGCCUCUGUGUGUUCCUGUAUGUGGGCGUGCACGUACCUGUCACAUGUGUACGCGCAGACCACAGGAUGUCCACACUGGCUUCCAAACACAUCU
    .(((.......(((((..((((...((((((.(.((.((((((((((.........)))))))))).)).)..))))))....))))..))))))))......... (-38.70)
    .(((...,,..(((((..((((.,.((((((.{.((.((((((((((.........)))))))))).)).}..))))))....))))..))))),,.......... [-40.43]
    ...........(((((..((((...((((((.(.((.((((((((((.........)))))))))).)).)..))))))....))))..)))))............ {-38.10 d=6.34}
    frequency of mfe structure in ensemble 0.0603635; ensemble diversity 9.71  
    # The above result extracts the first two lines, the id in the first line and the MFE value after the second line
    python 002getMFE.py 02RNAFold/mnv/mnv.res 02RNAFold/mnv/mnv_fix.res
    python 002getMFE.py 02RNAFold/snv/snv.res 02RNAFold/snv/snv_fix.res
    python 002getMFE.py 02RNAFold/ref/ref.res 02RNAFold/ref/ref_fix.res

# Calculate differences
    python 003getDiff_mnv.py 02RNAFold/ref/ref_fix.res 02RNAFold/mnv/mnv_fix.res 02RNAFold/mnv/mnv_diff.res 
    python 003getDiff_snv.py 02RNAFold/ref/ref_fix.res 02RNAFold/snv/snv_fix.res 02RNAFold/snv/snv_diff.res
    python 003getDiff_ms.py 02RNAFold/snv/snv_diff.res 02RNAFold/mnv/mnv_diff.res 02RNAFold/merge_res.txt