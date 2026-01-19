# 用途：批量把多套 SNV/MNV 数据转换为 BED/注释输入，调用 MNVAnno 与自编脚本生成基因/区域注释、区域计数、TSS 富集、UTR5/TSS 下游筛选等结果，供后续统计与绘图使用。
# 依赖：bedtools、Python(2/3) 环境以及 /home/caow/02mnv_new 与 /data/jinww/mnv/MNVAnno 下的脚本与数据库已可访问。
# 用法：直接执行即可；各输入路径已在脚本内写死，如需复用请调整数据路径与输出文件名。
# 注意：脚本包含大量一次性分析步骤（部分命令后台 & 执行），可按需求分段运行。
#   /data/jinww/mnv/analyse3/02data/03density/02datasets/snv_1000G.txt
# 在原始文件 raw.txt 的基础上生成 snv.tsv（TAB 分隔）
awk 'BEGIN{OFS="\t"}{
  split($1, v, ":");
  chr=v[1]; pos=v[2]; ref=v[3]; alt=v[4];
  id=$1; rs=".";
  print chr,pos,id,ref,alt,rs
}' /data/jinww/mnv/analyse3/02data/03density/02datasets/snv_UKB20w.txt > snv_UKB20w.anno.txt
bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i snv_1000G.anno.txt -t gene -D /data/jinww/mnv/MNVAnno/database/human/en_GRCh37.75_gene.anno.txt -o gene.txt &
bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i snv_UKB20w.anno.txt -t gene -D en_GRCh38.108_gene.anno.txt -o snv_UKB20w_gene.txt &

awk 'BEGIN{OFS="\t"}{
  split($1, v, ":");
  chr=v[1]; pos=v[2]; ref=v[3]; alt=v[4];
  id=$1; rs=".";
  print chr,pos,id,ref,alt,rs
}' /data/jinww/mnv/analyse3/02data/03density/02datasets/snv_GTEx.txt > snv_GTEx.anno.txt
python /home/caow/02mnv_new/scripts/SNVFineRegionAnno.py -i  ./snv_GTEx.anno.txt -D /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt   -o snv_GTEx.fine_region.tsv


# snv_TCGA.txt: rsID chr pos  ->  snv_TCGA.anno.txt: chr pos ID ref alt rsID
awk 'BEGIN{OFS="\t"}
{
  rs=$1; chr=$2; pos=$3;
  ref="N"; alt="N";
  print chr,pos,rs,ref,alt,rs
}' /data/jinww/mnv/analyse3/02data/03density/02datasets/snv_TCGA.txt > snv_TCGA.anno.txt
python /home/caow/02mnv_new/scripts/SNVFineRegionAnno.py -i  ./snv_TCGA.anno.txt -D /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt   -o snv_TCGA.fine_region.tsv

awk 'BEGIN{OFS="\t"}
{
  rs=$1; chr=$2; pos=$3;
  ref="N"; alt="N";
  print chr,pos,rs,ref,alt,rs
}' /data/jinww/mnv/analyse3/02data/03density/02datasets/snv_UKB50w.txt > snv_UKB50w.anno.txt
python /home/caow/02mnv_new/scripts/SNVFineRegionAnno.py -i  ./snv_UKB50w.anno.txt -D /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt   -o snv_UKB50w.fine_region.tsv


python /home/caow/02mnv_new/scripts/SNVFineRegionAnno.py -i  snv_UKB20w.anno.txt -D /data/jinww/mnv/MNVAnno/database/human/hg38_ATAC_aj.txt  -o snv_UKB20w.fine_atac.tsv


python3 /home/caow/02mnv_new/scripts/SNVBatchRegionCount.py \
  -i snv_TCGA.anno.txt \
  -o snv_TCGA.region.counts.tsv \
  --dbdir /data/jinww/mnv/MNVAnno/database/human \
  --build hg38 \
  --flank-to-intergenic

python3 /home/caow/02mnv_new/scripts/SNVBatchRegionCount.py \
  -i snv_UKB50w.anno.txt \
  -o snv_UKB50w.region.counts.tsv \
  --dbdir /data/jinww/mnv/MNVAnno/database/human \
  --build hg38 \
  --flank-to-intergenic

python3 /home/caow/02mnv_new/scripts/SNVBatchRegionCount.py \
  -i snv_UKB20w.anno.txt \
  -o snv_UKB20w.region.counts.tsv \
  --dbdir /data/jinww/mnv/MNVAnno/database/human \
  --build hg38 \
  --flank-to-intergenic

python3 /home/caow/02mnv_new/scripts/SNVBatchRegionCount.py \
  -i snv_1000G.anno.txt \
  -o snv_1000G.region.counts.tsv \
  --dbdir /data/jinww/mnv/MNVAnno/database/human \
  --build hg38 \
  --flank-to-intergenic

python3 /home/caow/02mnv_new/scripts/SNVBatchRegionCount.py \
  -i snv_GTEx.anno.txt \
  -o snv_GTEx.region.counts.tsv \
  --dbdir /data/jinww/mnv/MNVAnno/database/human \
  --build hg38 \
  --flank-to-intergenic

# 统计 SNV 和 MNV 在TSS区域的分布情况，来自https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3497849289_L2vPrabLMpdRIZnUAlnknOeKdlqp&clade=mammal&org=Human&db=hg38&hgta_group=regulation&hgta_track=robustPeaks&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C799%2C529-155%2C812%2C871&hgta_outputType=primaryTable&hgta_outFileName=
python3 /home/caow/02mnv_new/scripts/CountSNVMNVinRegions.py \
  -r ~/02mnv_new/addition/hgTables.txt \
  --snv snv_TCGA.anno.txt \
  --mnv /data/jinww/mnv/analyse3/02data/02datasets/TCGA.txt \
  -o tss_region.counts_TCGA.tsv

python3 /home/caow/02mnv_new/scripts/CountSNVMNVinRegions.py \
  -r ~/02mnv_new/addition/hgTables.txt \
  --snv snv_UKB50w.anno.txt \
  --mnv /data/jinww/mnv/analyse3/02data/02datasets/UKB50w.txt \
  -o tss_region.counts_UKB50w.tsv

python3 /home/caow/02mnv_new/scripts/CountSNVMNVinRegions.py \
  -r ~/02mnv_new/addition/hgTables.txt \
  --snv snv_UKB20w.anno.txt \
  --mnv /data/jinww/mnv/analyse3/02data/02datasets/UKB20w.txt \
  -o tss_region.counts_UKB20w.tsv

python3 /home/caow/02mnv_new/scripts/CountSNVMNVinRegions.py \
  -r ~/02mnv_new/addition/hgTables.txt \
  --snv snv_1000G.anno.txt \
  --mnv /data/jinww/mnv/analyse3/02data/02datasets/1000G.txt \
  -o tss_region.counts_1000G.tsv

python3 /home/caow/02mnv_new/scripts/CountSNVMNVinRegions.py \
  -r ~/02mnv_new/addition/hgTables.txt \
  --snv snv_GTEx.anno.txt \
  --mnv /data/jinww/mnv/analyse3/02data/02datasets/GTEx.txt \
  -o tss_region.counts_GTEx.tsv

# 修改了/home/caow/03mnv/analyse3/04adjustError/01region/01num中的数据统计，加入TSS数量统计


grep upstream_variant /home/caow/03mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline |grep 1000G |cut -f3 |sort -u >promotermnv1000G
grep upstream_variant /home/caow/03mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline |grep GTEx |cut -f3 |sort -u >promotermnvGTEx


awk -F'\t' '
BEGIN{OFS="\t"}
{
  chr=$1; poslist=$2; id=$3;
  n=split(poslist,a,",");
  for(i=1;i<=n;i++){
    p=a[i]+0;
    # BED: chr start end id
    # 单点用 [p-1, p)
    print "chr"chr, p-1, p, id;
  }
}' /home/caow/03mnv/analyse3/02data/02datasets/1000G.txt > 1000gmnv.pos.bed

awk -F'\t' 'BEGIN{OFS="\t"} {print $2,$3,$4,$5}' /home/caow/03mnv/analyse3/00addition/data/anno/hg38_TSS.txt > tss.bed

bedtools intersect -a 1000gmnv.pos.bed -b tss.bed -wa -wb \
| awk -F'\t' '
BEGIN{OFS="\t"}
{
  id=$4;
  tss=$8;   # intersect 输出：a的4列 + b的4列，所以b的label在第8列
  if(!( (id SUBSEP tss) in seen )){
    seen[id SUBSEP tss]=1;
    if(id in out) out[id]=out[id]","tss;
    else out[id]=tss;
  }
}
END{
  for(id in out) print id, out[id];
}' > 1000gMNV_to_TSS.tsv

awk -F'\t' '
BEGIN{OFS="\t"}
{
  chr=$1; poslist=$2; id=$3;
  n=split(poslist,a,",");
  for(i=1;i<=n;i++){
    p=a[i]+0;
    # BED: chr start end id
    # 单点用 [p-1, p)
    print "chr"chr, p-1, p, id;
  }
}' /home/caow/03mnv/analyse3/02data/02datasets/GTEx.txt > GTExmnv.pos.bed

awk -F'\t' 'BEGIN{OFS="\t"} {print $2,$3,$4,$5}' /home/caow/03mnv/analyse3/00addition/data/anno/hg38_TSS.txt > tss.bed

bedtools intersect -a GTExmnv.pos.bed -b tss.bed -wa -wb \
| awk -F'\t' '
BEGIN{OFS="\t"}
{
  id=$4;
  tss=$8;   # intersect 输出：a的4列 + b的4列，所以b的label在第8列
  if(!( (id SUBSEP tss) in seen )){
    seen[id SUBSEP tss]=1;
    if(id in out) out[id]=out[id]","tss;
    else out[id]=tss;
  }
}
END{
  for(id in out) print id, out[id];
}' > GTExMNV_to_TSS.tsv

# [1] "Baseline for 1000G : 0.041212660155921"
# [1] "Baseline for GTEx : 0.027862315966804"


cut -f3-5,7 /home/caow/03mnv/analyse3/02data/02datasets/1000G_mnv.txt >mnv1000g_pattern
cut -f3-5,7 /home/caow/03mnv/analyse3/02data/02datasets/GTEx_mnv.txt >mnvGTEx_pattern
cut -f1-4 /home/caow/03mnv/analyse3/03basePattern/01process/mnv.txt  >mnv1000gGTEx_pattern

grep upstream_variant /home/caow/03mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline |grep GTEx |cut -f3 |sort -u >promotermnvGTEx
grep "upstream_variant" /home/caow/03mnv/analyse3/02data/01origin/hg38mnv_addgene_oneline | \
grep -E "1000G|GTEx" | \
cut -f3 | sort -u > promotermnv_all

awk -F'\t' '
BEGIN{OFS="\t"}
{
  chr=$1; poslist=$2; id=$3;
  n=split(poslist,a,",");
  for(i=1;i<=n;i++){
    p=a[i]+0;
    # BED: chr start end id
    # 单点用 [p-1, p)
    print "chr"chr, p-1, p, id;
  }
}' /home/caow/03mnv/analyse3/02data/02datasets/1000G.txt > 1000gmnv.pos.bed

bedtools intersect -a GTExmnv.pos.bed -b tss.bed -wa -wb \
| awk -F'\t' '
BEGIN{OFS="\t"}
{
  id=$4;
  tss=$8;   # intersect 输出：a的4列 + b的4列，所以b的label在第8列
  if(!( (id SUBSEP tss) in seen )){
    seen[id SUBSEP tss]=1;
    if(id in out) out[id]=out[id]","tss;
    else out[id]=tss;
  }
}
END{
  for(id in out) print id, out[id];
}' > GTExMNV_to_TSS.tsv

# 1. 排序 (bedtools merge 要求输入必须按坐标排序)
# 2. 合并重叠区域
# 3. 计算 Start 和 End 的差值并求和
sort -k1,1 -k2,2n  promoter38.bed | \
bedtools merge -i stdin | \
awk '{sum += ($3 - $2)} END {print sum}'


# 1. 转换 MNV 文件为 BED 格式 (mnv.bed)
# 逻辑：跳过标题行，解析第6列(snvid)，找到最小和最大坐标
# 注意：BED 格式是 0-based，所以 Start = min_pos - 1
awk 'NR>1 {
    # 将第6列按逗号分割成多个 SNV
    n = split($6, snvs, ",");
    min_p = 1e15; max_p = 0;
    chrom = "";
    
    for (i=1; i<=n; i++) {
        # 再将每个 SNV 按冒号分割 (chr:pos:ref:alt)
        split(snvs[i], parts, ":");
        c = parts[1];
        p = parts[2];
        
        if (chrom == "") chrom = c; # 记录染色体
        if (p < min_p) min_p = p;   # 更新最小值
        if (p > max_p) max_p = p;   # 更新最大值
    }
    
    # 输出: Chr, Start-1, End, MNVID
    # 使用制表符分隔，确保 bedtools 兼容性
    printf("%s\t%d\t%d\t%s\n", chrom, min_p-1, max_p, $1);
}' /home/caow/03mnv/analyse3/03basePattern/01process/mnv.txt > mnv_1000ggtex.bed

# 2. 与 UTR5 文件取交集
# -wa: 输出原始的 A 文件内容 (即 MNV ID)
# -u: 只要有重叠就输出一次 (去重)
bedtools intersect -a mnv_1000ggtex.bed -b /home/caow/01miRNASNP/gene/utr5_cleaned.bed -wa | cut -f4 | sort -u > mnv_in_utr5_ids.txt

# 查看结果前几行
head mnv_in_utr5_ids.txt

bedtools intersect -a mnv_1000ggtex.bed -b /home/caow/03mnv/analyse3/00addition/data/tss_downstream_2000.bed -wa | cut -f4 | sort -u > mnv_in_tss_downstream_ids.txt

cut -f3 /home/caow/03mnv/analyse3/02data/02annotation/mature_miRNA.oneline.txt|sort -u >mirmnv
cut -f3 /home/caow/03mnv/analyse3/02data/02annotation/pre_miRNA.oneline.txt|sort -u >>mirmnv

for f in *.txt; do 
  awk '{print $2, $3-1, $4, $1}' OFS="\t" "$f" > "bed_output/${f%.txt}.bed"; 
done



awk -i inplace 'BEGIN{OFS="\t"} {$4="exon"; print}' bed_output/cds.bed
awk -i inplace 'BEGIN{OFS="\t"} {$4="UTR5"; print}' bed_output/utr5.bed
awk -i inplace 'BEGIN{OFS="\t"} {$4="UTR3"; print}' bed_output/utr3.bed
awk -i inplace 'BEGIN{OFS="\t"} {$4="intron"; print}' bed_output/intron.bed
awk -i inplace 'BEGIN{OFS="\t"} {$4="splice"; print}' bed_output/splice.bed
awk -i inplace 'BEGIN{OFS="\t"} {$4="intergenic"; print}' bed_output/intergenic.bed

awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1,$2,$3,$6}' /home/caow/03mnv/analyse3/00addition/data/constraint_z_genome_1kb.raw.download.txt > gnomad.expected.bed

# expected_share = expected_window * overlap_bp / window_len
# /home/caow/03mnv/analyse3/00addition/data
bedtools intersect -a gnomad.expected.bed -b tmp_anno_bed2 -wao \
| awk 'BEGIN{FS=OFS="\t"}
{
  a_chr=$1; a_s=$2; a_e=$3; expv=$4
  b_type=$8
  ov=$(NF)
  winlen=(a_e-a_s)
  if (winlen<=0) next
  if (ov<=0) next
  if (b_type=="." || b_type=="") next
  sum[b_type] += expv * (ov / winlen)
}
END{
  print "type","expected_sum"
  for (t in sum) printf "%s\t%.10f\n", t, sum[t]
}' | sort -k2,2nr > expected_by_type.hg38.tsv

awk '
{
    # 1. 跳过注释行
    if ($0 ~ /^#/ || $0 ~ /^track/ || $0 ~ /^browser/) next;

    # 2. 检查列数是否足够
    if (NF < 3) {
        print "Error: Line " NR " has fewer than 3 columns -> " $0;
        next;
    }

    # 3. 检查坐标是否为负数
    if ($2 < 0 || $3 < 0) {
        print "Error: Line " NR " contains negative coordinates -> " $0;
    }

    # 4. 检查 Start 是否小于 End (Start < End)
    if ($2 >= $3) {
        print "Error: Line " NR " has Start >= End -> " $0;
    }
}' tmp_anno_bed1

awk '($2 < $3 && $2 >= 0 && $3 >= 0)' tmp_anno_bed1 > tmp_anno_bed2


for f in $(ls /home/caow/03mnv/analyse3/00addition/data/bed_output/*.bed); do
  type=$(basename "$f" .bed)
  echo "Processing $type from $f"
  nohup python /home/caow/03mnv/analyse3/00addition/scripts/gnomad_snv_stat.py \
  --vcf /home/caow/03mnv/analyse3/00addition/data/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz \
  --bed $type=/home/caow/03mnv/analyse3/00addition/data/bed_output/$type.bed \
  --out /home/caow/03mnv/analyse3/00addition/results/gnomad_snv_counts_$type.tsv &
done


awk -F'\t' 'BEGIN{OFS="\t"}
$0!~/^#/ && $3=="gene"{
  chr=$1; start=$4-1; end=$5; strand=$7;

  id="."; name=".";

  # 提取 ID=... （如 gene:ENSG...）
  if ($9 ~ /(^|;)ID=/) {
    id=$9
    sub(/.*(^|;)ID=/, "", id)
    sub(/;.*/, "", id)
  }

  # 提取 Name=... （如 B3GALT6）
  if ($9 ~ /(^|;)Name=/) {
    name=$9
    sub(/.*(^|;)Name=/, "", name)
    sub(/;.*/, "", name)
  }

  print chr, start, end, id, name, strand
}' /home/caow/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3|awk 'BEGIN{OFS="\t"}{ if($1 !~ /^chr/) $1="chr"$1; print }'|sort -k1,1 -k2,2n > hg38genes.bed


sort -k1,1 -k2,2n /home/caow/03mnv/analyse3/00addition/data/1000gmnv.pos.bed > 1000gmnvpos.sorted.bed

bedtools intersect -a hg38genes.bed -b 1000gmnvpos.sorted.bed -wa -wb \
| awk 'BEGIN{OFS="\t"}
{
  gene=$4; name=$5; strand=$6; mnv=$10;   # 6列gene + 4列MNV => MNV_ID在第10列
  key = gene "|" name "|" strand "|" mnv
  gkey = gene "|" name "|" strand
  if (!(key in seen)) {
    seen[key]=1
    cnt[gkey]++
  }
}
END{
  for (g in cnt) {
    # 把 gkey 拆回三列输出
    n=split(g, a, "|")
    print a[1], a[2], a[3], cnt[g]
  }
}' \
| sort -k1,1 > 1000ggene_mnv_counts.tsv

awk 'BEGIN{OFS="\t"}
{
  id=$1; chr=$2; pos=$3;              # pos 是 1-based
  if (chr !~ /^chr/) chr="chr"chr;    # 统一 chr 前缀
  start=pos-1; end=pos;
  key=chr"\t"start"\t"end"\t"id;
  if(!(key in seen)){
    seen[key]=1;
    print chr, start, end, id
  }
}' /home/caow/03mnv/analyse3/02data/03density/02datasets/snv_1000G.txt | sort -k1,1 -k2,2n > 1000gsnv.sorted.bed

bedtools intersect -a hg38genes.bed -b 1000gsnv.sorted.bed -wa -wb \
| awk 'BEGIN{OFS="\t"}
{
  gene=$4; name=$5; strand=$6;
  snv=$10;                      # 6列gene + 4列snv => snv_id 在第10列
  key = gene "|" name "|" strand "|" snv
  gkey = gene "|" name "|" strand
  if (!(key in seen)) {
    seen[key]=1
    cnt[gkey]++
  }
}
END{
  for (g in cnt) {
    n=split(g,a,"|")
    print a[1], a[2], a[3], cnt[g]
  }
}' | sort -k1,1 > 1000ggene_snv_counts.tsv