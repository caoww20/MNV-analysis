## First step
# Perform GWAS analysis on cancer and extract GWAS significant loci (threshold 5e-6), annotate, calculate KM012 and plot
    sigGWA.5e-6.allcancer.results.xls
# Pick important points (survival curve can separate, KM significant)
    At least meet in one type can separate to get key.mnv
# Obtain rsid, and extract from GWAS
    cat key.mnv |while read id;do (grep $id /data/jinww/mnv/analyse3/02data/01origin/hg38mnv >>key.rsid);done &
    awk '{print $6}' key.rsid |sed 's/,/\n/g'|sort -u >key.rsid2
    cat key.rsid2 |while read id;do (echo $id && grep -w $id /data/jinww/04reference/publicDB/gwas/download/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv >>key.gwas);done &
# The genes where the above points fall, whether related to this cancer (literature search), manual
    
# Check if TCGA has this locus, whether the result is consistent

# Different from SNP, and significant loci in the region (regional effect map) 
    Temporarily not doable, because SNV data is not genotyped


# Merge all data
gwas:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/merge/gwas.txt
gene:/data/jinww/mnv/analyse3/02data/02annotation/gene_up.txt
gwas:/data/jinww/04reference/publicDB/gwas/download/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv
NCF:/data/jinww/04reference/publicDB/cancer_gene/NCG_cancerdrivers_annotation_supporting_evidence.tsv
intOGen:/data/jinww/04reference/publicDB/cancer_gene/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv
012:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/012/11Cancers.survival.p.txt
01:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/01/11Cancers.survival.p.txt
02:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/02/11Cancers.survival.p.txt
python merge.py
# Get real GWAS data
awk '$14 =="1" {print $3}' res.txt|sed 's/,/\n/g'|sort -u|while read id;do (echo $id && grep -w $id /data/jinww/04reference/publicDB/gwas/download/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv >>key.gwas);done &