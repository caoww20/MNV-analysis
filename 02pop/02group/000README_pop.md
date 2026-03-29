## Obtain Sample Information
    python 001getPopSample.py /data/jinww/04reference/publicDB/1000G/igsr-1000_genomes_30x_on_grch38.tsv /data/jinww/mnv/analyse3/02pop/02group/01origin/sample

## Capture MNVs from Shared SNVs
    # Extract MNVs from shared SNVs, result saved as MNVFromCommonSNV.txt
        python 001getMNVFromCommonSNV.py /data/jinww/mnv/analyse3/02pop/01overlap/venn/venn_1000G_common.txt /data/jinww/mnv/analyse3/02data/02datasets/1000G.txt 01origin/mnv.txt 02adjust/MNVFromCommonSNV.txt &

    # Obtain Genotype
        python /data/jinww/mnv/MNVAnno/tools/MNVGetGenotype.py -m MNVFromCommonSNV.txt -v /data/jinww/04reference/publicDB/1000G/hg38/IGSR.chr1.vcf -g 1 &

    # Convert to VCF format and adjust overlapping MNVs (i.e., +n)
        python 002toVCF.py 01origin/VCF.head 02adjust/MNVFromCommonSNV.genotype 02adjust/MNVFromCommonSNV.vcf &

## Capture Shared MNVs
    python 003getCommonMNV.py /data/jinww/mnv/analyse3/02pop/01overlap/venn/venn_1000G_common.txt 02adjust/MNVFromCommonSNV.vcf 02adjust/commonMNV.vcf &

## Calculate Allele Frequencies
    python 004getAF2.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv 02adjust/MNVFromCommonSNV.vcf 02adjust/MNVFromCommonSNV.freq
    python 004getAF2.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv 02adjust/commonMNV.vcf 02adjust/commonMNV.freq

## Perform Quality Control
    plink2 --vcf MNVFromCommonSNV.vcf --geno 0.05 --maf 0.01 --recode vcf --out MNVFromCommonSNVFilter &
    plink2 --vcf commonMNV.vcf --geno 0.05 --maf 0.01 --recode vcf --out commonMNVFilter &

## Calculate FST
    python 007getFSTSH.py MNVFromCommonSNVFilter.vcf AFR,AMR,EAS,EUR,SAS a.sh
    python 007getFSTSH.py 1000G.vcf AFR,AMR,EAS,EUR,SAS a.sh
    bash a.sh

    ## Merge Results
    python 008mergeFST.py ./03fst/ AFR,AMR,EAS,EUR,SAS 0 five_single_merge &  # "single" yields one result
    python 008mergeFST.py ./03fst/ AFR,AMR,EAS,EUR,SAS 1 five_region_merge &  # "region" yields two results: mean and weighted

## Count the number of MNVs with FST ≥ 0.15 in different populations (and those with FST > 0.15 in any population)
    Rscript 011getSigNum.R  # Generates sig.txt and population-specific significant data

## Use all.weir.fst to obtain SNV and MNV data with FST ≥ 0.15, extract genes containing at least 10 significant MNVs, and perform enrichment analysis
    awk '$3 >= 0.15 {print $1"\t"$2"\t"".\t""A\tT\t."}' 03fst/snv/all.weir.fst|grep -v CHROM|sed 's/chr//g' > 03fst/significant/snv_sig.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i snv_sig.txt -t gene -D en_GRCh38.108_gene.anno.txt -o snv_sig.gene

    awk '$3 >= 0.15 {print}' 03fst/filter/all_5.weir.fst > 03fst/significant/mnv_sig.txt &
    python 012getInfo2.py 03fst/significant/mnv_sig.txt 02adjust/MNVFromCommonSNVFilter.vcf /data/jinww/mnv/analyse2/02data/02annotation/gene.txt 03fst/significant/mnv_sig.gene &

    Rscript 013getGene.R