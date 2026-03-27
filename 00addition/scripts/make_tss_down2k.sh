#!/usr/bin/env bash
# Purpose: extract 2 kb downstream of transcript TSS from GFF3 and output BED6.
# Usage: ./make_tss_down2k.sh [GFF3] [OUT_BED]; default GFF3=$HOME/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3, output tss_downstream_2000.bed.
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

  # Parse ID or Name (avoid match third-arg for mawk compatibility)
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