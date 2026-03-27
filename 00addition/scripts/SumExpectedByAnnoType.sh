#!/usr/bin/env bash
# Purpose: sum gnomAD expected windows by BED annotation type (allocated by overlap proportion).
# Usage: ./SumExpectedByAnnoType.sh gnomad_windows.tsv[.gz] /path/to/anno_dir out.tsv
# Note: anno_dir should contain BED files; the script cats *.bed to aggregate.
set -euo pipefail

GNOMAD_TSV="${1:?Usage: $0 <gnomad_windows.tsv[.gz]> <anno_dir> <out.tsv>}"
ANNO_DIR="${2:?Usage: $0 <gnomad_windows.tsv[.gz]> <anno_dir> <out.tsv>}"
OUT_TSV="${3:?Usage: $0 <gnomad_windows.tsv[.gz]> <anno_dir> <out.tsv>}"

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }

tmp_gnomad_bed="$(mktemp)"
tmp_anno_bed="$(mktemp)"
trap 'rm -f "$tmp_gnomad_bed" "$tmp_anno_bed"' EXIT

# ---- 1) gnomAD: extract chr start end expected -> BED4
# Expected input column order: chrom start end element_id possible expected observed ...
# If your file has no header, remove NR==1{next}
# if [[ "$GNOMAD_TSV" == *.gz ]]; then
#   zcat "$GNOMAD_TSV"
# else
#   cat "$GNOMAD_TSV"
# fi | awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1,$2,$3,$6}' > "$tmp_gnomad_bed"

# ---- 2) Annotation: (type chrom start end) -> BED4: chr start end type
# Supports chrom values 1..22 or X/Y/MT or 23/24/25
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
#   # Expected format: type chrom start end ...
#   type=$1; chrom=$2; start=$3; end=$4
#   if (type=="" || chrom=="" || start=="" || end=="") next
#   chr=to_chr(chrom)
#   print chr,start,end,type
# }
# ' /dev/null > "$tmp_anno_bed"

# Merge common hg38 annotations in the directory (adjust as needed)
# Includes hg38_* and GRCh38 promoter (if they are type chrom start end format)
cat \
  "$ANNO_DIR"/*.bed \
  2>/dev/null  > "$tmp_anno_bed"

# ---- 3) Intersect and sum expected by type (allocated by overlap proportion)
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