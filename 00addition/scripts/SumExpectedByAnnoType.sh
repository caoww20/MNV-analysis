#!/usr/bin/env bash
# 功能：对 gnomAD expected 窗口按给定 BED 注释类型求期望总和（按重叠比例分摊）。
# 用法：./SumExpectedByAnnoType.sh gnomad_windows.tsv[.gz] /path/to/anno_dir out.tsv
# 说明：anno_dir 下需已准备好 BED，脚本会直接 cat *.bed 汇总后计算。
set -euo pipefail

GNOMAD_TSV="${1:?Usage: $0 <gnomad_windows.tsv[.gz]> <anno_dir> <out.tsv>}"
ANNO_DIR="${2:?Usage: $0 <gnomad_windows.tsv[.gz]> <anno_dir> <out.tsv>}"
OUT_TSV="${3:?Usage: $0 <gnomad_windows.tsv[.gz]> <anno_dir> <out.tsv>}"

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }

tmp_gnomad_bed="$(mktemp)"
tmp_anno_bed="$(mktemp)"
trap 'rm -f "$tmp_gnomad_bed" "$tmp_anno_bed"' EXIT

# ---- 1) gnomAD: 提取 chr start end expected -> BED4
# 期望输入列顺序类似：chrom start end element_id possible expected observed ...
# 若你的文件没有表头，把 NR==1{next} 去掉
# if [[ "$GNOMAD_TSV" == *.gz ]]; then
#   zcat "$GNOMAD_TSV"
# else
#   cat "$GNOMAD_TSV"
# fi | awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1,$2,$3,$6}' > "$tmp_gnomad_bed"

# ---- 2) 注释：把 (type chrom start end) -> BED4: chr start end type
# 支持 chrom 是 1..22 或 X/Y/MT 或 23/24/25
# awk -v OFS="\t" '
# function to_chr(c,   t){
#   t=c
#   if (t ~ /^chr/) return t
#   if (t=="23") t="X"
#   else if (t=="24") t="Y"
#   else if (t=="25") t="M"
#   else if (tolower(t)=="mt") t="M"
#   return "chr" t
# }
# BEGIN{FS=OFS="\t"}
# FNR==1{ }
# {
#   # 期望格式：type chrom start end ...
#   type=$1; chrom=$2; start=$3; end=$4
#   if (type=="" || chrom=="" || start=="" || end=="") next
#   chr=to_chr(chrom)
#   print chr,start,end,type
# }
# ' /dev/null > "$tmp_anno_bed"

# 把目录下常见 hg38 注释都并进来（按需增减）
# 这里包含 hg38_* 以及 GRCh38 promoter（如果它们也是 type chrom start end 格式）
cat \
  "$ANNO_DIR"/*.bed \
  2>/dev/null  > "$tmp_anno_bed"

# ---- 3) 相交并按 type 汇总 expected（按重叠比例分摊）
# expected_share = expected_window * overlap_bp / window_len
bedtools intersect -a "$tmp_gnomad_bed" -b "$tmp_anno_bed" -wao \
| awk 'BEGIN{FS=OFS="\t"}
{
  a_chr=$1; a_s=$2; a_e=$3; exp=$4
  b_type=$8
  ov=$(NF)
  winlen=(a_e-a_s)
  if (winlen<=0) next
  if (ov<=0) next
  if (b_type=="." || b_type=="") next
  sum[b_type] += exp * (ov / winlen)
}
END{
  print "type","expected_sum"
  for (t in sum) printf "%s\t%.10f\n", t, sum[t]
}' | sort -k2,2nr > "$OUT_TSV"

echo "Wrote: $OUT_TSV" >&2