# ==Calculate MNV for each method=================================================
# MNVAnno (base)
    # Get MNV
    bash 00cost.sh "bash /data/jinww/mnv/MNVAnno/MNVIdentify.sh -i all_res2.vcf" "VCFCheck|MNVIdentify"

# Qing (MNVAnno2)
    bash 00cost.sh "python get_mnv2.py all_res2_pass.vcf" "hail"

    # Convert to the following format
    #chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
    python 01fixFormat.py all_res2_pass.mnv qing.mnv

# BCFtools (base)
    bash 00cost.sh "bcftools csq -f /data/jinww/04reference/genome/Human/ensemble38/Homo_sapiens.GRCh38.dna.toplevel.fa -g Homo_sapiens.GRCh38.109.chr.gff3 -o res.txt -O v --threads 1  --ncsq 63 all_res2.vcf" bcftools
    # --ncsq Maximum number of consequences per haplotype per site, can be adjusted as needed
    # MNV results need further extraction, and AC values are not calculated, more focused on annotation
    # Also can only get MNVs in CDS regions @ symbol !!! Can identify spaced MNVs and 3-base MNVs

    # Convert to the following format
    #chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
    grep '@' res.txt > temp
    python 01fixFormat.py res.txt temp bcftools.mnv 
    rm temp

# MACARON (py27)
    ls *markdup.vcf|while read id;do (bash 00cost.sh "bash 00MACARON.sh $id" "gatk|snpEff|MACARON_fix2");done
    # Input files do not need to be sorted
    # Can only process one VCF at a time, only identifies MNVs causing missense, can have gaps
    # Due to this reason, maximum memory usage does not increase with sample number
    # Also can only get MNVs in CDS regions @ symbol !!! Can identify spaced MNVs and 3-base MNVs

    # Convert to the following format
    #chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
    python 01fixFormat.py file temp.mnv # This code is simply used for current results, not applicable to new data
    ls *txt|while read id;do (python 01fixFormat.py  $id temp.mnv);done 
    sort -u temp.mnv > MACARON.mnv
    rm temp.mnv

# vardict (base)
    <!-- Documentation: https://github.com/AstraZeneca-NGS/VarDictJava -->
    awk -F '\t' '$1 == 22 && $3 == "gene" {OFS="\t"; print $1, $4, $5}'  /data/jinww/mnv/MNVAnno/database/human_build/ensembl/hg38/Homo_sapiens.GRCh38.108.chr.gtf >gene_chr22.bed
    ls /data/jinww/mnvcompare/data/*_markdup.bam|while read id;do (bash 00cost.sh "perl vardict -G /data/jinww/04reference/genome/Human/hg38.fa -b $id -f 0.01 -c 1 -S 2 -E 3  gene_chr22.bed | teststrandbias.R | var2vcf_valid.pl -E -f 0.01 > ${id%_*}_vardict.vcf &" "vardict|teststrandbias|var2vcf_valid");done
    # -g 4 can be omitted
    # Because there may be errors during running, starting with some bed regions having no data, so use for loop as follows
    # One sample at a time
    for ((i=0; i<4039; i++)); do (sed -n "$((i+1)),$((i+1))p" my.bed > temp.bed;echo $i; perl vardict -G /data/jinww/04reference/genome/Human/hg38.fa -b /data/jinww/mnvcompare/data/HG00524_markdup.bam -f 0.01 -c 1 -S 2 -E 3 temp.bed | ../teststrandbias.R | ../var2vcf_valid.pl -E -f 0.01 >> temp.vcf); done
    
    for ((i=0; i<4039; i++)); do (sed -n "$((i+1)),$((i+1))p" my.bed > ./HG00663/temp.bed;echo $i; perl vardict -G /data/jinww/04reference/genome/Human/hg38.fa -b /data/jinww/mnvcompare/data/HG00663_markdup.bam -f 0.01 -c 1 -S 2 -E 3 ./HG00663/temp.bed | teststrandbias.R | var2vcf_valid.pl -E -f 0.01 >> ./HG00663/temp.vcf); done &
    # Record time and memory
    bash 00test.sh &
    ls HG*/temp.vcf|while read id;do (grep -v "#" $id| grep Complex|grep PASS > ${id%/*}/${id%/*}.mnv);done
    ls NA*/temp.vcf|while read id;do (grep -v "#" $id| grep Complex|grep PASS > ${id%/*}/${id%/*}.mnv);done
    # Can only identify consecutive MNVs !!! Confirmed, tested with aaa_2.vcf
    # Convert to the following format
    #chr	pos	MNVID	ref	alt	rsID	distance	MNVType	AC	adjust_AC	adjust_AF
    ls HG*/*.mnv|while read id;do (python 01fixFormat.py $id  ${id%.*}_fix.mnv);done
    ls NA*/*.mnv|while read id;do (python 01fixFormat.py $id  ${id%.*}_fix.mnv);done
    cat HG*/*_fix.mnv >a
    cat NA*/*_fix.mnv >b
    cat a b >vardict.mnv;rm a b
    awk -F'\t' '{sub(/^[^\t]+\t/, ""); print}' vardict.mnv|sort -u >vardict_simple.mnv
    
# MAC (root environment)
    # Extract sites from VCF
    # grep -v '#' ../../data/HG00524_markdup.vcf|awk 'length($4) <= 1 && length($5) <= 1 {print $1, $2, $4, $5}'|sed 's/ /./g'
    ls ../../data/*_markdup.vcf|while read id;do (grep -v '#' $id|awk 'length($4) <= 1 && length($5) <= 1 {print $1, $2, $4, $5}'|sed 's/ /./g' > /data/jinww/mnvcompare/method/MAC/${id##*/} &);done
    # Run MAC_v1.0.pl
    <!-- perl MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -->
    perl MAC_v1.2.pl -i ./data/HG00524_markdup.vcf -bam ../../data/HG00524_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa
    ls *vcf |while read id;do ( bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i $id -bam /data/jinww/mnvcompare/data/${id%.*}.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl");done
    # Considering it runs very slowly, submit 10 at a time
    # Convert to the following format
    ls *_counts.txt|while read id;do (python 01fixFormat.py $id  MAC.mnv );done

# Run MAC separately for each sample, due to time cost, each bam file takes 276m
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG00524_markdup.vcf -bam /data/jinww/mnvcompare/data/HG00524_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG00663_markdup.vcf -bam /data/jinww/mnvcompare/data/HG00663_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG00693_markdup.vcf -bam /data/jinww/mnvcompare/data/HG00693_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG01046_markdup.vcf -bam /data/jinww/mnvcompare/data/HG01046_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG01286_markdup.vcf -bam /data/jinww/mnvcompare/data/HG01286_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG01444_markdup.vcf -bam /data/jinww/mnvcompare/data/HG01444_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG01817_markdup.vcf -bam /data/jinww/mnvcompare/data/HG01817_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02253_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02253_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02265_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02265_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02277_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02277_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02922_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02922_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02924_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02924_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02946_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02946_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG02948_markdup.vcf -bam /data/jinww/mnvcompare/data/HG02948_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG03073_markdup.vcf -bam /data/jinww/mnvcompare/data/HG03073_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG03616_markdup.vcf -bam /data/jinww/mnvcompare/data/HG03616_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG03617_markdup.vcf -bam /data/jinww/mnvcompare/data/HG03617_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG03642_markdup.vcf -bam /data/jinww/mnvcompare/data/HG03642_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG03706_markdup.vcf -bam /data/jinww/mnvcompare/data/HG03706_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG04235_markdup.vcf -bam /data/jinww/mnvcompare/data/HG04235_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i NA12286_markdup.vcf -bam /data/jinww/mnvcompare/data/NA12286_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i NA12829_markdup.vcf -bam /data/jinww/mnvcompare/data/NA12829_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i NA12874_markdup.vcf -bam /data/jinww/mnvcompare/data/NA12874_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i NA20517_markdup.vcf -bam /data/jinww/mnvcompare/data/NA20517_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i NA20800_markdup.vcf -bam /data/jinww/mnvcompare/data/NA20800_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &

# ==Analysis=================================================
# valid
    example: samtools view -hb -L sub1.bed sample02.bam > sub1.chr22_21349676-21349677.sample02.bam
    sub1.bed
    chr22	21349676
    The naming format of output BAM file should be the same.
    Once subset BAM files are generated, run MACARON_validate.sh.
    example: MACARON_validate.sh sub1.chr22_21349676-21349677.sample02.bam
    # 改写成自己的脚本 
    cat sample.list |while read id;do (bash 01my_validate.sh /data/jinww/mnvcompare/data/ $id &);done &
    # 考虑到运算时间，拆分成5份跑
    cat sample.list |while read id;do (bash 01my_validate.sh /data/jinww/mnvcompare/data/ $id part1.txt &);done &
    cat sample.list |while read id;do (cat part*/${id}_read.txt > merge/${id}_merge.txt);done
# Extract MNVs present in VCF from total
    


# View sam or bam flag
https://www.omicshare.com/tools/Home/Soft/samflag
https://broadinstitute.github.io/picard/explain-flags.html
The seq in bam files are all consistent with ref, see /data/jinww/mnvcompare/data/test_sam/