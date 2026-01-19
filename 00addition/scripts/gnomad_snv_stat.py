"""
功能：统计 gnomAD VCF 在多个 BED 区域内的 SNP 数量与覆盖长度，输出汇总 TSV。
用法示例：
    python gnomad_snv_stat.py --vcf gnomad.vcf.bgz --bed UTR5=utr5.bed --bed exon=exon.bed --out snv_counts.tsv
或：
    python gnomad_snv_stat.py --vcf gnomad.vcf.bgz --bed-list bed_list.txt --out snv_counts.tsv
依赖：pysam；VCF 需有索引(.tbi/.csi)。
"""

import pysam
import os
import sys
import argparse

class VariantCounter:
    def __init__(self, vcf_path):
        """
        初始化：加载 VCF 文件（只需加载一次）
        """
        print(f"正在加载 VCF 索引: {vcf_path} ...")
        self.vcf = pysam.VariantFile(vcf_path, "r")
        self.vcf_path = vcf_path

    def _merge_intervals(self, intervals):
        """
        内部辅助函数：合并重叠区间
        输入: list of [start, end]
        输出: list of [start, end] (已合并重叠部分)
        """
        if not intervals:
            return []

        # 按起始位置排序
        intervals.sort(key=lambda x: x[0])

        merged = []
        for current in intervals:
            if not merged:
                merged.append(current)
            else:
                last = merged[-1]
                # 如果当前区间的开始 <= 上一个区间的结束，说明重叠
                if current[0] < last[1]:
                    # 合并：结束位置取两者的最大值
                    last[1] = max(last[1], current[1])
                else:
                    merged.append(current)
        return merged

    def _read_and_merge_bed(self, bed_path):
        """
        读取 BED 文件，按染色体分组，并合并每一组的重叠区间
        """
        chrom_intervals = {}  # 结构: {'chr1': [[100,200], [300,400]], 'chr2': ...}

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
                    continue  # 跳过格式错误的行

        # 对每个染色体的区间进行合并
        for chrom in chrom_intervals:
            chrom_intervals[chrom] = self._merge_intervals(chrom_intervals[chrom])

        return chrom_intervals

    def count_snps_in_region_type(self, bed_path, region_name="Unknown"):
        """
        核心对外接口：统计指定 BED 文件覆盖区域内的 SNP 总数
        返回: (total_snps, total_bp_covered)
        """
        print(f"正在处理区域类型: [{region_name}] (文件: {bed_path})...")

        # 1. 读取并合并区间，消除重叠
        merged_intervals_map = self._read_and_merge_bed(bed_path)

        total_snps = 0
        total_bp_covered = 0

        # 2. 遍历每个染色体的合并区间去查询 VCF
        for chrom, intervals in merged_intervals_map.items():
            for start, end in intervals:
                total_bp_covered += (end - start)
                try:
                    for record in self.vcf.fetch(chrom, start, end):
                        if len(record.ref) == 1 and record.alts:
                            if any(len(alt) == 1 for alt in record.alts):
                                total_snps += 1
                except ValueError:
                    # 染色体不存在或命名风格不匹配
                    pass

        print(f" -> 结果: [{region_name}] 覆盖长度 {total_bp_covered} bp, 包含 {total_snps} 个 SNP")
        return total_snps, total_bp_covered


def _parse_beds_from_args(bed_args, bed_list_path):
    """
    支持两种输入：
    1) --bed NAME=PATH 可重复
    2) --bed-list 两列文件：name<TAB>path（也接受空格分隔）
    返回 dict: {name: path}
    """
    bed_files = {}

    if bed_args:
        for item in bed_args:
            if "=" not in item:
                raise ValueError(f"--bed 参数格式应为 NAME=PATH，收到: {item}")
            name, path = item.split("=", 1)
            name = name.strip()
            path = path.strip()
            if not name or not path:
                raise ValueError(f"--bed 参数格式应为 NAME=PATH，收到: {item}")
            bed_files[name] = path

    if bed_list_path:
        with open(bed_list_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"--bed-list 每行至少两列 name path，收到: {line}")
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
        description="统计 gnomAD VCF 在多个 BED 区域中的 SNP 数量，并输出汇总表"
    )
    parser.add_argument(
        "--vcf",
        required=False,
        default="/home/caow/03mnv/analyse3/00addition/data/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
        help="输入 VCF(.vcf.gz/.bgz) 路径（需存在索引 .tbi/.csi）",
    )
    parser.add_argument(
        "--bed",
        action="append",
        default=[],
        help="BED 定义（可重复）：NAME=PATH，例如 --bed Enhancers=/path/enhancer.bed",
    )
    parser.add_argument(
        "--bed-list",
        default=None,
        help="BED 列表文件：每行 name<tab>path（或空格分隔），# 开头为注释",
    )
    parser.add_argument(
        "--out",
        default="snv_counts.tsv",
        help="输出汇总 TSV 文件路径（默认 snv_counts.tsv）",
    )

    args = parser.parse_args()

    vcf_file = args.vcf
    if not os.path.exists(vcf_file):
        print("错误: 找不到 VCF 文件，请检查路径。", file=sys.stderr)
        sys.exit(1)

    try:
        bed_files = _parse_beds_from_args(args.bed, args.bed_list)
    except Exception as e:
        print(f"错误: 解析 BED 参数失败：{e}", file=sys.stderr)
        sys.exit(2)

    if not bed_files:
        print("错误: 未提供任何 BED。请使用 --bed 或 --bed-list。", file=sys.stderr)
        sys.exit(2)

    counter = VariantCounter(vcf_file)

    results_rows = []
    for region_name, bed_path in bed_files.items():
        if not os.path.exists(bed_path):
            print(f"错误: 找不到文件 {bed_path}", file=sys.stderr)
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
        print("错误: 没有任何可写出的结果（BED 都不存在或处理失败）。", file=sys.stderr)
        sys.exit(3)

    _write_results_tsv(args.out, results_rows)

    print("\n最终统计汇总:")
    for r in results_rows:
        print(f"{r['region']}: {r['snps']} SNPs, {r['bp_covered']} bp, density={r['snps_per_bp']}")
    print(f"\n已写出结果到: {args.out}")