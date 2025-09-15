## Error correction =========================
## Count several indicators
    # Number of datasets
        gnomAD      5977564     # grep GnomAD /data/jinww/mnv/analyse3/02data/01origin/hg38mnv|wc -l 
        FiveData    3984258     # grep -E 'TCGA|GTEx|1000G|UKB50w|UKB20w' /data/jinww/mnv/analyse3/02data/01origin/hg38mnv|wc -l 
        total       8199654     # wc -l /data/jinww/mnv/analyse3/02data/01origin/hg38mnv 
    # Number of MNVs falling in CDS
        gnomAD      68174   # grep exon /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene|grep GnomAD|wc -l
        FiveData    52200   # grep exon /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene|grep -E 'TCGA|GTEx|1000G|UKB50w|UKB20w'|wc -l
        total       102483   # grep exon /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene| wc -l
# All transcripts
    # Number of MNVs falling in amino acids
        gnomAD      33889   # grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene|grep GnomAD|wc -l
        FiveData    10730    # grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene|grep -E 'TCGA|GTEx|1000G|UKB50w|UKB20w'|wc -l
        total       39487   # grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene|wc -l     
    # Number of MNVs falling in amino acids and changing amino acids completely different from SNVs (only count 2~3 joint MNVs within the same codon)
        # Extract data
            grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline|grep GnomAD > 01allTrans/gnomAD.txt &
            grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline|grep -E 'TCGA|GTEx|1000G|UKB50w|UKB20w' > 01allTrans/fiveData.txt &
            grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline > 01allTrans/total.txt &
        # Extract only 2~3 joint MNVs within the same codon
            python 001codonMNV.py gnomAD.txt gnomAD_fix.mnv &
            python 001codonMNV.py fiveData.txt fiveData_fix.mnv &
            python 001codonMNV.py total.txt total_fix.mnv &
        # Number of MNVs falling in amino acids (2~3 joints)
            total   37207   awk '{print $3}' total_fix.mnv |sort -u |wc -l
        # Extract SNV data
            python 002codonSNV.py gnomAD_fix.mnv gnomAD_fix.snv &
            python 002codonSNV.py fiveData_fix.mnv fiveData_fix.snv &
            python 002codonSNV.py total_fix.mnv total_fix.snv &
        # Annotate SNVs
            bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i gnomAD_fix.snv -t gene -D en_GRCh38.108_gene.anno.txt -o gnomAD_fix_snv_anno.txt &
            bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i fiveData_fix.snv -t gene -D en_GRCh38.108_gene.anno.txt -o fiveData_fix_snv_anno.txt &
            bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i total_fix.snv -t gene -D en_GRCh38.108_gene.anno.txt -o total_fix_snv_anno.txt &
            python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i gnomAD_fix_snv_anno.txt -f T &
            python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i fiveData_fix_snv_anno.txt -f T &
            python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i total_fix_snv_anno.txt -f T &
        ## Add SNV results and base changes
            python 003annotateSNV.py gnomAD_fix_snv_anno.oneline.txt /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt gnomAD_fix.snv.new &
            python 003annotateSNV.py fiveData_fix_snv_anno.oneline.txt /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt fiveData_fix.snv.new &
            python 003annotateSNV.py total_fix_snv_anno.oneline.txt /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt total_fix.snv.new &
        ## Integrate SNV and MNV
            python 004merge.py gnomAD_fix.snv.new  gnomAD_fix.mnv merge_gnomAD.txt &
            python 004merge.py fiveData_fix.snv.new  fiveData_fix.mnv merge_fiveData.txt &
            python 004merge.py total_fix.snv.new  total_fix.mnv merge_total.txt &
        ## Calculate the number of MNVs falling in amino acids and changing amino acids completely different from SNVs (as long as one transcript's MNV amino acid change is completely different from SNV)
            Rscript 005summary.R
            [1] "gnomAD"
                [1] "mnvid: 17663"
                [1] "gene: 9274"
                [1] "transcript: 27790"
            [1] "fiveData"
                [1] "mnvid: 4378"
                [1] "gene: 3214"
                [1] "transcript: 8485"
            [1] "total"
                [1] "mnvid: 19388"
                [1] "gene: 9726"
                [1] "transcript: 28964"
        ## Calculate how many SNVs are misannotated
            python 006getSNVMisAnnoted.py merge_total.txt
            53981 70872 76.21%

# Canonical transcripts
    # Number of MNVs falling in amino acids
        gnomAD      32070   # grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_canonical|grep GnomAD|wc -l
        FiveData    9481    # grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_canonical|grep -E 'TCGA|GTEx|1000G|UKB50w|UKB20w'|wc -l
        total       36761   # grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_canonical|wc -l
    # Number of MNVs falling in amino acids and changing amino acids completely different from SNVs (only count 2~3 joint MNVs within the same codon)
        # Extract data
            grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_canonical_oneline|grep GnomAD > 02canonicalTrans/gnomAD.txt &
            grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_canonical_oneline|grep -E 'TCGA|GTEx|1000G|UKB50w|UKB20w' > 02canonicalTrans/fiveData.txt &
            grep "@" /data/jinww/mnv/analyse3/02data/01origin/hg38mnv_addgene_canonical_oneline > 02canonicalTrans/total.txt &
        # Extract only 2~3 joint MNVs within the same codon
            python 001codonMNV.py gnomAD.txt gnomAD_fix.mnv &
            python 001codonMNV.py fiveData.txt fiveData_fix.mnv &
            python 001codonMNV.py total.txt total_fix.mnv &
        # Extract SNV data
            python 002codonSNV.py gnomAD_fix.mnv gnomAD_fix.snv &
            python 002codonSNV.py fiveData_fix.mnv fiveData_fix.snv &
            python 002codonSNV.py total_fix.mnv total_fix.snv &
        # Annotate SNVs
            bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i gnomAD_fix.snv -t gene -D en_GRCh38.108_gene.anno.txt -o gnomAD_fix_snv_anno.txt &
            bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i fiveData_fix.snv -t gene -D en_GRCh38.108_gene.anno.txt -o fiveData_fix_snv_anno.txt &
            bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i total_fix.snv -t gene -D en_GRCh38.108_gene.anno.txt -o total_fix_snv_anno.txt &
            python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i gnomAD_fix_snv_anno.txt -f T &
            python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i fiveData_fix_snv_anno.txt -f T &
            python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i total_fix_snv_anno.txt -f T &
        ## Add SNV results and base changes
            python 003annotateSNV.py gnomAD_fix_snv_anno.oneline.txt /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt gnomAD_fix.snv.new &
            python 003annotateSNV.py fiveData_fix_snv_anno.oneline.txt /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt fiveData_fix.snv.new &
            python 003annotateSNV.py total_fix_snv_anno.oneline.txt /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt total_fix.snv.new &
        ## Integrate SNV and MNV
            python 004merge.py gnomAD_fix.snv.new  gnomAD_fix.mnv merge_gnomAD.txt &
            python 004merge.py fiveData_fix.snv.new  fiveData_fix.mnv merge_fiveData.txt &
            python 004merge.py total_fix.snv.new  total_fix.mnv merge_total.txt &
        ## Calculate the number of MNVs falling in amino acids and changing amino acids completely different from SNVs (as long as one transcript's MNV amino acid change is completely different from SNV)
            Rscript 005summary.R
            [1] "gnomAD"
                [1] "mnvid: 16684"
                [1] "gene: 8760"
                [1] "transcript: 8787"
            [1] "fiveData"
                [1] "mnvid: 3946"
                [1] "gene: 2905"
                [1] "transcript: 2915"
            [1] "total"
                [1] "mnvid: 18155"
                [1] "gene: 9156"
                [1] "transcript: 9185"
    # Flowchart drawing
        python 005alluvial.py merge_total.txt alluvial.txt
## Get the number of each gene causing amino acid MNVs
    python 007getAA.py /data/jinww/mnv/analyse2/02data/02annotation/gene.oneline.txt ./
    awk '{print $1}' AA_perGene.txt |sort -u|wc -l # 36790 MNVs can cause amino acid changes
    awk '{print $1}' start_change_perGene.txt |sort -u|wc -l # 2002 MNVs can cause start codon changes, 249 cause start_lost 
    awk '{print $1}' stop_change_perGene.txt |sort -u|wc -l # 2510 MNVs can cause stop codon changes


# There is an important note here
    # 2	171324260,171324261	MNV02434689	C,A	T,G	rs372978109,.	ENSG00000123600 METTL8 ENST00000375258 exon10,exon10 G1136A,T1135C @*379Q stop_lost High
    # 2	171324260	.	C	T	.	ENSG00000123600 METTL8 ENST00000375258 exon10 G1136A R379Q missense_variant Moderate
    # 2	171324261	.	A	G	.	ENSG00000123600 METTL8 ENST00000375258 exon10 T1135C *379R stop_lost High
    # Theoretically all should be *, but after checking the data, the reference base at position 61 is G, not A, so during MNV identification, it was done using the real ref, so this is different from the single point results, these can be ignored.
