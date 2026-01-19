## 00addition/scripts 说明

本目录包含补充分析脚本，主要围绕 MNV/SNV 区域注释、期望值计算、TFBS/miRNA 下游分析及绘图。路径默认指向本机数据目录，复用时请按需修改。

### 环境依赖
- `bash`、`bedtools`（sort/merge/intersect/subtract）
- Python 3：`pysam`（解析 VCF），常规标准库
- R：`tidyverse`/`dplyr`、`ggplot2`、`reshape2`、`gridExtra`、`stringr`、`GenomicRanges`、`rtracklayer`、`ggpubr`

### 脚本概览
- `workflow.sh`：长管线，批量把多套 SNV/MNV 数据转 BED/注释，调用 MNVAnno 与自编脚本生成基因/区域注释、区域计数、TSS 富集、UTR5/TSS 下游筛选等结果。直接运行即可，必要时修改脚本内硬编码路径与输出名。
- `anno2bed.py`：将当前目录下 `hg38_*.txt`（列：type chr start end …）转换为 BED4 (`chr start-1 end type`)，输出到 `data/bed_output/`。示例：`python anno2bed.py`。
- `gnomad_snv_stat.py`：统计 gnomAD VCF 在多个 BED 区域的 SNP 数量与覆盖长度。示例：
	`python gnomad_snv_stat.py --vcf gnomad.vcf.bgz --bed UTR5=utr5.bed --bed exon=exon.bed --out snv_counts.tsv`
- `GnomadExpectedByGeneRegion.py`：按 SNVBatchRegionCount 优先级把 gnomAD expected 窗口分摊到 splice/exon/UTR/intron/up/down/intergenic。示例：
	`python GnomadExpectedByGeneRegion.py -i gnomad.windows.tsv.gz -o expected_by_region.tsv --gene-db en_GRCh38.108_canonical_gene.anno.txt --up 2000 --down 1000 --splice 2`
- `make_tss_down2k.sh`：从 GFF3 提取转录本 TSS 下游 2kb，输出 BED6。示例：`./make_tss_down2k.sh Homo_sapiens.GRCh38.110.gff3 tss_downstream_2000.bed`。
- `SumExpectedByAnnoType.sh`：对 gnomAD expected 窗口按给定 BED 注释类型求期望总和（按重叠比例分摊）。示例：`./SumExpectedByAnnoType.sh gnomad_windows.tsv anno_dir expected_by_type.tsv`。
- `plot_bak.R`：汇总 MNV/SNV 区域分布、期望比率及多项下游分析的绘图脚本，需先在脚本内指定的工作目录准备好 `01num`、`expected_by_type` 等输入。
- `tfbs.r`：基于已加载的 `df`、`expected_var`、`merged_data` 等数据，计算期望/长度与 GC、重复序列比例并绘图。
- `tfbs1.r`：分析 MNV- TFBS 结果与靶基因约束力（LOEUF/pLI），需要 `tfbs_res`、`jaspar_meta`、gnomAD 约束表及 GFF3（默认 `~/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3`）。
- `xinjibing12snvenrich.R`：评估心肌病 GWAS 信号在高 FST MNV 区域的富集，使用泊松/二项/超几何检验。

### 运行提示
- 脚本内大量路径为绝对路径（`/home/caow`, `/data/jinww/mnv` 等），迁移到新环境请先统一数据位置。
- `workflow.sh` 含多段一次性分析与后台命令，可按需求逐段运行并确认输入已就绪。
- 需要 R 脚本时，确保所需数据框已预先读入，或按脚本注释调整读取路径。
