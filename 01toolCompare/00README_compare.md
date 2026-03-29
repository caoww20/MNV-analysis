# ==Calculate MNV for each method=================================================
# MNVAnno (base)
    # Get MNV
    bash 00cost.sh "bash /data/jinww/mnv/MNVAnno/MNVIdentify.sh -i all_res2.vcf" "VCFCheck|MNVIdentify"

# Qing (MNVAnno2)
    bash 00cost.sh "python get_mnv2.py all_res2_pass.vcf" "hail"

# BCFtools (base)
    bash 00cost.sh "bcftools csq -f /data/jinww/04reference/genome/Human/ensemble38/Homo_sapiens.GRCh38.dna.toplevel.fa -g Homo_sapiens.GRCh38.109.chr.gff3 -o res.txt -O v --threads 1  --ncsq 63 all_res2.vcf" bcftools

# MACARON (py27)
    ls *markdup.vcf|while read id;do (bash 00cost.sh "bash 00MACARON.sh $id" "gatk|snpEff|MACARON_fix2");done

# vardict (base)
    <!-- Documentation: https://github.com/AstraZeneca-NGS/VarDictJava -->    
    for ((i=0; i<4039; i++)); do (sed -n "$((i+1)),$((i+1))p" my.bed > ./HG00663/temp.bed;echo $i; perl vardict -G /data/jinww/04reference/genome/Human/hg38.fa -b /data/jinww/mnvcompare/data/HG00663_markdup.bam -f 0.01 -c 1 -S 2 -E 3 ./HG00663/temp.bed | teststrandbias.R | var2vcf_valid.pl -E -f 0.01 >> ./HG00663/temp.vcf); done &
    
    
# MAC (root environment)
    bash /data/jinww/mnvcompare/method/MAC/00cost.sh "perl /data/jinww/mnvcompare/method/MAC/MAC_v1.2.pl -i HG00524_markdup.vcf -bam /data/jinww/mnvcompare/data/HG00524_markdup.bam -r /data/jinww/04reference/genome/Human/hg38.fa &" "MAC_v1.2.pl" &
