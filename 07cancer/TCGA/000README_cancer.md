## TCGA different cancers MNV counts
    grep TCGA /data/jinww/mnv/analyse3/02data/01origin/hg19mnv > TCGA.mnv
## Annotate results
    # Since it is version 19, annotate separately here
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i TCGA.mnv -t gene -D en_GRCh37.75_gene.anno.txt -o gene.txt &
## TCGA eQTL MNVs
    # Download data from pancanMNVeQTL to get eQTL.txt
    # Convert eQTL.txt results to MNVID results
    python 001transID.py TCGA.txt gene.txt eQTL.txt eQTL_fix.txt
## How many cancers can each eQTLMNV affect, most eQTL-MNVs have tissue specificity, but a few can affect multiple cancers
    # Use R's table function to count each MNV's cancer number based on rsid
## Draw circo plot, different regions MNVs affecting cancer numbers (promotor, UTR5, CDS, splice, UTR3, intro, other)
    R
## Draw different regions beta values (including trans?)
    R
## Overlap of these eQTLMNV genes with cancer genes
    python 003getGene.py eQTL_fix.txt eQTL_fix.gene
    /data/jinww/04reference/publicDB/cancer_gene/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv
    /data/jinww/04reference/publicDB/cancer_gene/NCG_cancerdrivers_annotation_supporting_evidence.tsv
## Extract GWAS results overlap with these known loci and LD close MNVs
    # Only focus on lung
    # Overlap lung eQTLMNV with lung GWAS results, see falling in promotor, UTR5, CDS, splice, UTR3
    python 002getSig.py eQTL_fix.txt UKB50wMNVGWAS/03.Result/01.allCancer/14.lung/fastGWA_GLMM_final.fastGWA eQTL_gwas.txt
## Example 【Different from SNP, and effective！！】
