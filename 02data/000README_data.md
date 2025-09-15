## Preprocessing (copy data from Luo Haohui, this version retains MHC and black)
    ## MHC
    https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
    https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
    hg19 chr6:28477797-33448354
    hg38 chr6:28510120-33480577
## Copy data from Luo Haohui
    cp /home/luohh/MNV_analysis/Statistic3/all_dataset/hg19/all.hg19.results.mnv ./
    cp /home/luohh/MNV_analysis/Statistic3/all_dataset/hg38/all.hg38.results.mnv ./
## Preprocess data (hg38 version processing)
    awk '{if ($8<6) print $o}' all.hg19.results.mnv >hg19mnv &
    awk '{if ($8<6) print $o}' all.hg38.results.mnv >hg38mnv &
## Split datasets into individual datasets and mark MNV sources
    python 001splitData.py 01origin/hg38mnv 02datasets
## Annotate hg38_humanMNV
    # Gene annotation
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t gene -D en_GRCh38.108_gene.anno.txt -o gene.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t gene -D en_GRCh38.108_canonical_gene.anno.txt -o gene_canonical.txt &
    # Non-coding annotation
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t non -D hg38_circRNA.txt -o circ.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t non -D hg38_lncRNA.txt -o lnc.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t non -D hg38_pre_miRNA.txt -o pre_miRNA.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t non -D hg38_mature_miRNA.txt -o mature_miRNA.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t non -D hg38_piRNA.txt -o piR.txt &
    # Other library region annotation
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D hg38_ATAC_aj.txt -o atac.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D hg38_CE_aj.txt -o CE.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D hg38_enhancer_aj.txt -o enhancer.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D hg38_miRBS_aj.txt -o miRBS.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i hg38mnv -t reg -D hg38_TFBS_aj.txt -o TFBS.txt &
    # Get oneline results
    ls *txt |while read id;do (python /data/jinww/mnv/MNVAnno/tools/toOneLine.py -i $id &);done
# Later need to merge data with population and annotation information, so use scripts to merge here
    python 002mergePopAnno.py 01origin/hg38mnv 02annotation/gene.txt 0 01origin/hg38mnv_addgene &
    python 002mergePopAnno.py 01origin/hg38mnv 02annotation/gene.oneline.txt 1 01origin/hg38mnv_addgene_oneline &
    python 002mergePopAnno.py 01origin/hg38mnv 02annotation/gene_canonical.txt 0 01origin/hg38mnv_addgene_canonical &
    python 002mergePopAnno.py 01origin/hg38mnv 02annotation/gene_canonical.oneline.txt 1 01origin/hg38mnv_addgene_canonical_oneline &
## Calculate the relationship between SNV density and MNV density
    SNV: Use 01getSNV.py
        python 01getSNV.py /data/jinww/04reference/publicDB/1000G/hg38/ 1000G snv_1000G.txt
        python 01getSNV.py /data/jinww/04reference/publicDB/GTEx/hg38/ GTEx snv_GTEx.txt
        python 01getSNV.py /data/jinww/04reference/publicDB/PancanQTL/vcf/ TCGA snv_TCGA.txt
        sort -u snv_TCGA.txt >a;rm snv_TCGA.txt;mv a snv_TCGA.txt
        python 01getSNV.py /data/jinww/04reference/publicDB/UKB20w/hg38/ UKB20w snv_UKB20.txt &
        python 01getSNV.py /data/jinww/04reference/publicDB/UKB50w/ UKB50w snv_UKB50.txt &
        # For snv_UKB50, to speed up, use 01getSNV_forUKB50w.py to parallel all chromosomes splitChr, then merge results
    MNV: Use 02getMNV.py
        python 001getMNV.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv 1000G mnv_1000G.txt &
        python 001getMNV.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv GTEx mnv_GTEx.txt &
        python 001getMNV.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv UKB20w mnv_UKB20w.txt &
        python 001getMNV.py /data/jinww/mnv/analyse3/02data/01origin/hg19mnv TCGA mnv_TCGA.txt &
        python 001getMNV.py /data/jinww/mnv/analyse3/02data/01origin/hg19mnv UKB50w mnv_UKB50w.txt &
    ## Calculate density, TCGA chr18, UKB20w chr3
        python 003getDes.py 02datasets/snv_1000G.txt hg38 03density/snv_1000G_Density.txt &
        python 003getDes.py 02datasets/snv_GTEx.txt hg38 03density/snv_GTEx_Density.txt &
        python 003getDes.py 02datasets/snv_UKB20w.txt hg38 03density/snv_UKB20w_Density.txt &
        python 003getDes.py 02datasets/snv_TCGA.txt hg19 03density/snv_TCGA_Density.txt &
        python 003getDes.py 02datasets/snv_UKB50w.txt hg19 03density/snv_UKB50w_Density.txt &

        python 003getDes.py 02datasets/mnv_1000G.txt hg38 03density/mnv_1000G_Density.txt &
        python 003getDes.py 02datasets/mnv_GTEx.txt hg38 03density/mnv_GTEx_Density.txt &
        python 003getDes.py 02datasets/mnv_UKB20w.txt hg38 03density/mnv_UKB20w_Density.txt &
        python 003getDes.py 02datasets/mnv_TCGA.txt hg19 03density/mnv_TCGA_Density.txt &
        python 003getDes.py 02datasets/mnv_UKB50w.txt hg19 03density/mnv_UKB50w_Density.txt &
## Why choose MNVs with 5 joints or less
    Data comes from statistic.file5, this part is mainly handled by Luo Haohui
## Count the number of each dataset and their mutual overlap
    plot.R
    ## Statistics of identification results ## Count the total number of MNVs and descriptive information
## Count how many MNV AC errors were corrected in total, how many MNVs were deleted that were mistakenly identified as low-joint MNVs
    awk '{print $1,$2,$4,$5,$9,$10}' /home/luohh/UKB50wMNVIdentify/hg19/UKB50w.all.mnv > hg19_UKB50w.all.mnv &
    awk '{print $1,$2,$4,$5,$9,$10}' /home/luohh/UKB20wMNVIdentify/hg38/UKB20w.all.mnv >hg38_UKB20w.all.mnv &
    ls /home/luohh/TCGAmnv_identify/hg19/*_phase/Result/*.all.mnv |while read id;do (awk '{print $1,$2,$4,$5,$9,$10}' $id > hg19_${id##*/});done &
    awk '{print $1,$2,$4,$5,$9,$10}' /home/luohh/GTExmnv_identify/hg38/GTEx_838/MNV_identify/GTEx.all.mnv >hg38_GTEx.all.mnv &
    awk '{print $1,$2,$4,$5,$9,$10}' /home/luohh/1000Kmnv_identify/hg38/1000K_3202/1000G.all.mnv > hg38_1000G.all.mnv &
    # Total number
    ls *mnv >mnv.list
    cat mnv.list |while read id;do (wc -l $id |cut -d' ' -f1 >> total.txt );done
    # How many were corrected
    cat mnv.list |while read id;do (awk '$5!=$6 {print $o}' $id |wc -l >> adjust.txt);done
    # How many were deleted
    cat mnv.list |while read id;do (awk '$6=="0" {print $o}' $id |wc -l >>filter.txt);done
    # Merge
    paste mnv.list total.txt adjust.txt filter.txt > merge.txt

## R program to view data for each dataset
    df<-fread('flag_datasets',header = T,sep='\t',stringsAsFactors = F,data.table = F)
    df$gnomad<-df$GnomAD_WGS+df$GnomAD_WES
    df[df$gnomad==2,'gnomad']<-1
    df$diff<-df$Sum-df$gnomad
    df<-df[df$diff>0,]
    nrow(df)
    nrow(df[df$mnvtype==2,])
    nrow(df[df$mnvtype==3,])
    nrow(df[df$mnvtype==4,])
    nrow(df[df$mnvtype==5,])