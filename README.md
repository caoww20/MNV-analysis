# MNV-analysis

This is a script-first analysis repository for multiple nucleotide variants (MNVs). The project is organized as staged pipelines, covering data preparation, population differences, base-level patterns, functional impact, ncRNA, TFBS, and cancer-related analyses.

## 1. Project Scope

- Focus: distribution, functional differences, and biological impact of MNVs and SNVs.
- Organization: numbered stages (00-07), each with standalone scripts and docs.
- Stack: Python + R + Bash.
- Characteristics: data-driven workflows, multiple external dependencies, no single package entrypoint, and no unified test suite.

## 2. Repository Structure

```text
MNV-analysis/
├─ 00addition/      # Supplementary annotation, region statistics, expected values, and plotting
├─ 01toolCompare/   # MNV tool comparison
├─ 02data/          # Dataset split/merge, source tagging, SNV/MNV density calculation
├─ 02pop/           # Population overlap/grouping, FST, and enrichment
├─ 03basePattern/   # Base-level pattern processing
├─ 04adjustError/   # Regional and functional error adjustment
├─ 05ncRNA/         # pre-miRNA, lncRNA, and miRNA binding-site workflows
├─ 06TFBS/          # Promoter and TFBS impact analyses
├─ 07cancer/        # Cancer-related analyses
└─ LICENSE
```

## 3. Quick Start

### 3.1 Execution Principles

There is no single global command for this repository. Run workflows from each module directory according to its README.

Common command patterns:

```bash
python <script>.py <args...>
Rscript <script>.R
bash <script>.sh
```

### 3.2 Recommended Reading Order

1. [00addition/README.md](00addition/README.md)
2. [02data/000README_data.md](02data/000README_data.md)
3. [02pop/02group/000README_pop.md](02pop/02group/000README_pop.md)
4. [03basePattern/000README_basePattern.md](03basePattern/000README_basePattern.md)
5. [04adjustError/01region/000README_adjustError.md](04adjustError/01region/000README_adjustError.md)
6. [04adjustError/02function/000README_adjustError.md](04adjustError/02function/000README_adjustError.md)
7. [05ncRNA/01pre_miRNA/000README_pre_miRNA.md](05ncRNA/01pre_miRNA/000README_pre_miRNA.md)
8. [05ncRNA/02lncRNA/000README_lncRNA.md](05ncRNA/02lncRNA/000README_lncRNA.md)
9. [06TFBS/000README_TFBS.md](06TFBS/000README_TFBS.md)
10. [01toolCompare/00README_compare.md](01toolCompare/00README_compare.md)

## 4. Module Navigation

### 00addition

A collection of supplementary scripts: annotation conversion, gnomAD statistics, expected-value calculation, region/TSS analyses, and plotting.

- Entry doc: [00addition/README.md](00addition/README.md)
- Typical entry script: [00addition/scripts/workflow.sh](00addition/scripts/workflow.sh)

### 01toolCompare

MNV tool comparison, including performance and output benchmarking.

- Doc: [01toolCompare/00README_compare.md](01toolCompare/00README_compare.md)

### 02data

Raw data split, annotation merge, and SNV/MNV density calculation; this is the core upstream stage for downstream analyses.

- Doc: [02data/000README_data.md](02data/000README_data.md)

### 02pop

Includes overlap analysis (01overlap) and grouped population analysis (02group), covering VCF processing, AF/FST calculation, and enrichment.

- Overlap doc: [02pop/01overlap/000README_pop.md](02pop/01overlap/000README_pop.md)
- Group doc: [02pop/02group/000README_pop.md](02pop/02group/000README_pop.md)

### 03basePattern

Base-level SNV/MNV feature processing and pattern analysis.

- Doc: [03basePattern/000README_basePattern.md](03basePattern/000README_basePattern.md)

### 04adjustError

Contains both regional and functional parts, focusing on MNV/SNV comparison, error adjustment, and coding-impact evaluation.

- Regional doc: [04adjustError/01region/000README_adjustError.md](04adjustError/01region/000README_adjustError.md)
- Functional doc: [04adjustError/02function/000README_adjustError.md](04adjustError/02function/000README_adjustError.md)

### 05ncRNA

Includes workflows for pre-miRNA, lncRNA, and miRNA binding sites (UTR3).

- pre-miRNA doc: [05ncRNA/01pre_miRNA/000README_pre_miRNA.md](05ncRNA/01pre_miRNA/000README_pre_miRNA.md)
- lncRNA doc: [05ncRNA/02lncRNA/000README_lncRNA.md](05ncRNA/02lncRNA/000README_lncRNA.md)
- miRBS doc: [05ncRNA/03miRBS_UTR3/000README_miRBS_UTR3.md](05ncRNA/03miRBS_UTR3/000README_miRBS_UTR3.md)
- RNAsnp supplementary doc: [05ncRNA/02lncRNA/02RNAsnp/README_RNAsnp.md](05ncRNA/02lncRNA/02RNAsnp/README_RNAsnp.md)

### 06TFBS

Promoter extraction, TFBS prediction, gain/loss statistics, and downstream matching analyses.

- Doc: [06TFBS/000README_TFBS.md](06TFBS/000README_TFBS.md)

### 07cancer

Cancer-related analysis module (including TCGA-related scripts and subdirectories).

## 5. Environment and Dependencies

Common dependencies (used by different modules):

- Python 3 (some legacy scripts may mix versions)
- R (common packages: ggplot2, dplyr, reshape2, stringr, GenomicRanges, rtracklayer, etc.)
- Bash
- bedtools
- plink2
- bcftools (in selected workflows)
- MNVAnno (core annotation tool)
- ViennaRNA/RNAfold (used in some ncRNA workflows)

Linux is recommended. If you run on Windows, using WSL with unified path mapping is strongly recommended.

## 6. License

See [LICENSE](LICENSE).
