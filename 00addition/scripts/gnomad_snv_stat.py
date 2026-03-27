"""
Purpose: count SNPs and covered length for gnomAD VCF across multiple BED regions and output a summary TSV.
Example usage:
    python gnomad_snv_stat.py --vcf gnomad.vcf.bgz --bed UTR5=utr5.bed --bed exon=exon.bed --out snv_counts.tsv
Or:
    python gnomad_snv_stat.py --vcf gnomad.vcf.bgz --bed-list bed_list.txt --out snv_counts.tsv
Dependency: pysam; VCF must have an index (.tbi/.csi).
"""

import pysam
import os
import sys
import argparse

class VariantCounter:
    def __init__(self, vcf_path):
        """
        Initialize: load VCF once
        """
        print(f"Loading VCF index: {vcf_path} ...")
        self.vcf = pysam.VariantFile(vcf_path, "r")
        self.vcf_path = vcf_path

    def _merge_intervals(self, intervals):
        """
        Internal helper to merge overlapping intervals
        Input: list of [start, end]
        Output: list of [start, end] (merged)
        """
        if not intervals:
            return []

        # Sort by start position
        intervals.sort(key=lambda x: x[0])

        merged = []
        for current in intervals:
            if not merged:
                merged.append(current)
            else:
                last = merged[-1]
                # If current start <= previous end, intervals overlap
                if current[0] < last[1]:
                    # Merge: take max end
                    last[1] = max(last[1], current[1])
                else:
                    merged.append(current)
        return merged

    def _read_and_merge_bed(self, bed_path):
        """
        Read BED, group by chromosome, and merge overlaps per group
        """
        chrom_intervals = {}  # Structure: {'chr1': [[100,200], [300,400]], 'chr2': ...}

        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(('#', 'track', 'browser')):
                    continue

                parts = line.split()
                chrom = parts[0]
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                try:
                    start = int(parts[1])
                    end = int(parts[2])

                    if chrom not in chrom_intervals:
                        chrom_intervals[chrom] = []
                    chrom_intervals[chrom].append([start, end])
                except ValueError:
                    continue  # Skip malformed lines

        # Merge intervals per chromosome
        for chrom in chrom_intervals:
            chrom_intervals[chrom] = self._merge_intervals(chrom_intervals[chrom])

        return chrom_intervals

    def count_snps_in_region_type(self, bed_path, region_name="Unknown"):
        """
        Main API: count total SNPs covered by the given BED regions
        Returns: (total_snps, total_bp_covered)
        """
        print(f"Processing region type: [{region_name}] (file: {bed_path})...")

        # 1. Read and merge intervals to remove overlaps
        merged_intervals_map = self._read_and_merge_bed(bed_path)

        total_snps = 0
        total_bp_covered = 0

        # 2. Query VCF using merged intervals per chromosome
        for chrom, intervals in merged_intervals_map.items():
            for start, end in intervals:
                total_bp_covered += (end - start)
                try:
                    for record in self.vcf.fetch(chrom, start, end):
                        if len(record.ref) == 1 and record.alts:
                            if any(len(alt) == 1 for alt in record.alts):
                                total_snps += 1
                except ValueError:
                    # Chromosome missing or naming mismatch
                    pass

        print(f" -> Result: [{region_name}] covered length {total_bp_covered} bp, contains {total_snps} SNPs")
        return total_snps, total_bp_covered


def _parse_beds_from_args(bed_args, bed_list_path):
    """
    Supports two input modes:
    1) --bed NAME=PATH (repeatable)
    2) --bed-list two-column file: name<TAB>path (space-delimited also accepted)
    Returns dict: {name: path}
    """
    bed_files = {}

    if bed_args:
        for item in bed_args:
            if "=" not in item:
                raise ValueError(f"--bed format should be NAME=PATH, got: {item}")
            name, path = item.split("=", 1)
            name = name.strip()
            path = path.strip()
            if not name or not path:
                raise ValueError(f"--bed format should be NAME=PATH, got: {item}")
            bed_files[name] = path

    if bed_list_path:
        with open(bed_list_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"--bed-list requires at least two columns (name path), got: {line}")
                name = parts[0]
                path = parts[1]
                bed_files[name] = path

    return bed_files


def _write_results_tsv(out_path, results_rows):
    """
    results_rows: list of dict with keys:
      region, bed_path, bp_covered, snps, snps_per_bp
    """
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as w:
        w.write("region\tbed_path\tbp_covered\tsnps\tsnps_per_bp\n")
        for r in results_rows:
            w.write(
                f"{r['region']}\t{r['bed_path']}\t{r['bp_covered']}\t{r['snps']}\t{r['snps_per_bp']}\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count SNPs in gnomAD VCF across multiple BED regions and output a summary table"
    )
    parser.add_argument(
        "--vcf",
        required=False,
        default="/home/caow/03mnv/analyse3/00addition/data/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
        help="Input VCF (.vcf.gz/.bgz) path (index .tbi/.csi required)",
    )
    parser.add_argument(
        "--bed",
        action="append",
        default=[],
        help="BED definition (repeatable): NAME=PATH, e.g. --bed Enhancers=/path/enhancer.bed",
    )
    parser.add_argument(
        "--bed-list",
        default=None,
        help="BED list file: each line name<tab>path (or space-delimited); # starts comments",
    )
    parser.add_argument(
        "--out",
        default="snv_counts.tsv",
        help="Output summary TSV path (default snv_counts.tsv)",
    )

    args = parser.parse_args()

    vcf_file = args.vcf
    if not os.path.exists(vcf_file):
        print("Error: VCF file not found. Check the path.", file=sys.stderr)
        sys.exit(1)

    try:
        bed_files = _parse_beds_from_args(args.bed, args.bed_list)
    except Exception as e:
        print(f"Error: failed to parse BED arguments: {e}", file=sys.stderr)
        sys.exit(2)

    if not bed_files:
        print("Error: no BED provided. Use --bed or --bed-list.", file=sys.stderr)
        sys.exit(2)

    counter = VariantCounter(vcf_file)

    results_rows = []
    for region_name, bed_path in bed_files.items():
        if not os.path.exists(bed_path):
            print(f"Error: file not found {bed_path}", file=sys.stderr)
            continue

        snps, bp = counter.count_snps_in_region_type(bed_path, region_name=region_name)
        snps_per_bp = (snps / bp) if bp else 0.0
        results_rows.append(
            {
                "region": region_name,
                "bed_path": bed_path,
                "bp_covered": bp,
                "snps": snps,
                "snps_per_bp": f"{snps_per_bp:.6e}",
            }
        )

    if not results_rows:
        print("Error: no results to write (BEDs missing or failed).", file=sys.stderr)
        sys.exit(3)

    _write_results_tsv(args.out, results_rows)

    print("\nFinal summary:")
    for r in results_rows:
        print(f"{r['region']}: {r['snps']} SNPs, {r['bp_covered']} bp, density={r['snps_per_bp']}")
    print(f"\nWrote results to: {args.out}")