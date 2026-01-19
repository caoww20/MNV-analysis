#!/usr/bin/env bash
# 功能：从 GFF3 提取转录本 TSS 下游 2kb 区域，输出 BED6。
# 用法：./make_tss_down2k.sh [GFF3] [OUT_BED] ；默认 GFF3=$HOME/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3，输出 tss_downstream_2000.bed。
set -euo pipefail

GFF3="${1:-$HOME/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3}"
OUT_BED="${2:-tss_downstream_2000.bed}"

awk -F'\t' -v OFS='\t' '
BEGIN { n=2000 }
$0 ~ /^#/ { next }
($3=="transcript" || $3=="mRNA") {
  chr=$1
  s1=$4+0
  e1=$5+0
  strand=$7
  attr=$9

  # 解析 ID 或 Name（不用 match 的第三参数，兼容 mawk）
  name="."
  if (attr ~ /ID=/) {
    name=attr
    sub(/^.*ID=/, "", name)
    sub(/;.*$/, "", name)
  } else if (attr ~ /Name=/) {
    name=attr
    sub(/^.*Name=/, "", name)
    sub(/;.*$/, "", name)
  }

  if (strand=="+") {
    start0 = s1 - 1
    end   = start0 + n
  } else if (strand=="-") {
    tss   = e1
    start0 = tss - n
    end    = tss
  } else next

  if (start0 < 0) start0 = 0
  if (end < 0) next

  print chr, start0, end, name, 0, strand
}
' "$GFF3" > "$OUT_BED"

echo "Wrote: $OUT_BED" >&2