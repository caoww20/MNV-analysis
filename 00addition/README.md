## 00addition/scripts Notes

This directory contains supplementary analysis scripts, mainly for MNV/SNV region annotation, expected-value calculation, TFBS/miRNA downstream analysis, and plotting. Paths default to local data directories; update them as needed when reusing.

### Environment dependencies
- `bash`, `bedtools` (sort/merge/intersect/subtract)
- Python 3: `pysam` (for VCF parsing), standard library
- R: `tidyverse`/`dplyr`, `ggplot2`, `reshape2`, `gridExtra`, `stringr`, `GenomicRanges`, `rtracklayer`, `ggpubr`

### Script overview
- `workflow.sh`: A long pipeline that batch-converts multiple SNV/MNV datasets to BED/annotations, calls MNVAnno and custom scripts to generate gene/region annotations, region counts, TSS enrichment, UTR5/TSS downstream filters, and related outputs. Run directly; edit hard-coded paths and output names if needed.
- `anno2bed.py`: Converts `hg38_*.txt` in the current directory (columns: type chr start end ...) to BED4 (`chr start-1 end type`) and writes to `data/bed_output/`. Example: `python anno2bed.py`.
- `gnomad_snv_stat.py`: Counts SNPs and coverage length for gnomAD VCF across multiple BED regions. Example:
	`python gnomad_snv_stat.py --vcf gnomad.vcf.bgz --bed UTR5=utr5.bed --bed exon=exon.bed --out snv_counts.tsv`
- `GnomadExpectedByGeneRegion.py`: Allocates gnomAD expected windows by SNVBatchRegionCount priority into splice/exon/UTR/intron/up/down/intergenic. Example:
	`python GnomadExpectedByGeneRegion.py -i gnomad.windows.tsv.gz -o expected_by_region.tsv --gene-db en_GRCh38.108_canonical_gene.anno.txt --up 2000 --down 1000 --splice 2`
- `make_tss_down2k.sh`: Extracts transcript TSS downstream 2 kb from GFF3 and outputs BED6. Example: `./make_tss_down2k.sh Homo_sapiens.GRCh38.110.gff3 tss_downstream_2000.bed`.
- `SumExpectedByAnnoType.sh`: Sums expected values for gnomAD windows by BED annotation type (allocated by overlap proportion). Example: `./SumExpectedByAnnoType.sh gnomad_windows.tsv anno_dir expected_by_type.tsv`.
- `plot_bak.R`: Plotting script for MNV/SNV region distributions, expected ratios, and multiple downstream analyses; prepare inputs like `01num`, `expected_by_type` in the working directory defined in the script.
- `tfbs.r`: Based on loaded data like `df`, `expected_var`, `merged_data`, computes expectation/length and GC/repeat proportions, then plots.
- `tfbs1.r`: Analyzes MNV-TFBS results and target-gene constraint (LOEUF/pLI); requires `tfbs_res`, `jaspar_meta`, gnomAD constraint tables, and GFF3 (default `~/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3`).
- `xinjibing12snvenrich.R`: Evaluates enrichment of cardiomyopathy GWAS signals in high-FST MNV regions using Poisson/binomial/hypergeometric tests.

### Running notes
- Many paths in scripts are absolute (e.g., `/home/caow`, `/data/jinww/mnv`); unify data locations before moving to a new environment.
- `workflow.sh` includes multiple one-off analyses and background commands; run section by section as needed and verify inputs are ready.
- For R scripts, ensure required data frames are preloaded or adjust read paths per script comments.
