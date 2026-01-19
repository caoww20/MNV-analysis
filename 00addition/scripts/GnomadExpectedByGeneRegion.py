#!/usr/bin/env python3
"""
功能：按与 SNVBatchRegionCount 相同的优先级，将 gnomAD expected 窗口按基因区域类型（splice/exon/UTR/intron/up/down/intergenic）分摊汇总。
用法示例：
    python GnomadExpectedByGeneRegion.py -i gnomad.windows.tsv.gz -o expected_by_region.tsv \
        --gene-db /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_canonical_gene.anno.txt \
        --up 2000 --down 1000 --splice 2 --coord bed0 --expected-col 5
依赖：bedtools 可用；输入列需包含 chrom start end expected(默认第6列)。
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# ------------------------
# 与 SNVBatchRegionCount.py 一致的 transcript 解析（1-based inclusive）
# ------------------------
@dataclass(frozen=True)
class Transcript:
    tid: str
    chr: str
    strand: str
    tx_start: int
    tx_end: int
    cds_start: int
    cds_end: int
    exons: Tuple[Tuple[int, int], ...]

def norm_chr(c: str) -> str:
    c = c.strip()
    return c[3:] if c.startswith("chr") else c

def to_chr(c: str) -> str:
    c = norm_chr(c)
    if c in ("MT", "Mt", "mt"):
        c = "M"
    return "chr" + c

def is_int(s: str) -> bool:
    try:
        int(s)
        return True
    except Exception:
        return False

def parse_int_list(s: str) -> List[int]:
    if s is None:
        return []
    s = s.strip().rstrip(",")
    if not s:
        return []
    sep = "," if "," in s else (";" if ";" in s else None)
    if sep is None:
        return [int(s)] if is_int(s) else []
    out = []
    for x in s.split(sep):
        x = x.strip()
        if x and is_int(x):
            out.append(int(x))
    return out

def parse_exons(exon_count: int, exon_starts_field: str, exon_ends_field: str) -> Optional[Tuple[Tuple[int, int], ...]]:
    if exon_count == 1 and is_int(exon_starts_field) and is_int(exon_ends_field):
        s = int(exon_starts_field); e = int(exon_ends_field)
        if s > e:
            s, e = e, s
        return ((s, e),)

    starts = parse_int_list(exon_starts_field)
    ends = parse_int_list(exon_ends_field)
    if not starts or not ends or len(starts) != len(ends):
        return None
    exons = []
    for s, e in zip(starts, ends):
        if s > e:
            s, e = e, s
        exons.append((s, e))
    exons.sort()
    return tuple(exons)

def load_transcripts(gene_db_path: str) -> Dict[str, List[Transcript]]:
    per_chr: Dict[str, List[Transcript]] = defaultdict(list)
    with open(gene_db_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 12:
                continue
            if not (is_int(a[5]) and is_int(a[6]) and is_int(a[7]) and is_int(a[8]) and is_int(a[9])):
                continue

            tid = a[0]
            chrid = to_chr(a[1])
            strand = a[2] if a[2] in ("+", "-") else "+"
            tx_start = int(a[5]); tx_end = int(a[6])
            cds_start = int(a[7]); cds_end = int(a[8])
            exon_count = int(a[9])

            exons = parse_exons(exon_count, a[10], a[11])
            if exons is None:
                continue

            if tx_start > tx_end:
                tx_start, tx_end = tx_end, tx_start
            if cds_start > cds_end:
                cds_start, cds_end = cds_end, cds_start

            per_chr[chrid].append(Transcript(
                tid=tid, chr=chrid, strand=strand,
                tx_start=tx_start, tx_end=tx_end,
                cds_start=cds_start, cds_end=cds_end,
                exons=exons
            ))
    for c in per_chr:
        per_chr[c].sort(key=lambda t: (t.tx_start, t.tx_end, t.tid))
    return per_chr

# ------------------------
# 坐标工具：1-based inclusive -> BED(0-based, half-open)
# ------------------------
def to_bed0(s1: int, e1: int) -> Tuple[int, int]:
    # 输入 1-based inclusive [s1,e1] -> 输出 bed [s0,e0)
    return (s1 - 1, e1)

def clip_1based(s: int, e: int, lo: int, hi: int) -> Optional[Tuple[int, int]]:
    if s > e:
        s, e = e, s
    s = max(lo, s)
    e = min(hi, e)
    if s <= e:
        return (s, e)
    return None

def intersect_1based(a: Tuple[int, int], b: Tuple[int, int]) -> Optional[Tuple[int, int]]:
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    if s <= e:
        return (s, e)
    return None

# ------------------------
# 生成各类 region 的原始 interval（可能重叠），后续用 bedtools 按优先级做互斥化
# ------------------------
def emit_raw_regions(per_chr: Dict[str, List[Transcript]], up: int, down: int, splice: int, outdir: str) -> Dict[str, str]:
    paths = {k: os.path.join(outdir, f"raw.{k}.bed") for k in ["splice","exon","UTR5","UTR3","intron","upstream","downstream"]}
    fps = {k: open(p, "w", encoding="utf-8") for k,p in paths.items()}

    def write_iv(kind: str, chrid: str, s1: int, e1: int):
        s0, e0 = to_bed0(s1, e1)
        if s0 < 0:
            s0 = 0
        if e0 > s0:
            fps[kind].write(f"{chrid}\t{s0}\t{e0}\n")

    for chrid, ts in per_chr.items():
        for t in ts:
            tx_lo, tx_hi = t.tx_start, t.tx_end

            # upstream/downstream（按链，均为 1-based inclusive）
            if t.strand == "+":
                up_iv = clip_1based(t.tx_start - up, t.tx_start - 1, 1, 10**12)
                down_iv = clip_1based(t.tx_end + 1, t.tx_end + down, 1, 10**12)
            else:
                up_iv = clip_1based(t.tx_end + 1, t.tx_end + up, 1, 10**12)
                down_iv = clip_1based(t.tx_start - down, t.tx_start - 1, 1, 10**12)
            if up_iv:
                write_iv("upstream", chrid, up_iv[0], up_iv[1])
            if down_iv:
                write_iv("downstream", chrid, down_iv[0], down_iv[1])

            # UTR 区域（基于 CDS 边界，落在 exon 内）
            # 先定义 genomic 范围（仍是 1-based inclusive）
            if t.strand == "+":
                utr5_rng = (tx_lo, t.cds_start - 1)
                utr3_rng = (t.cds_end + 1, tx_hi)
            else:
                utr3_rng = (tx_lo, t.cds_start - 1)
                utr5_rng = (t.cds_end + 1, tx_hi)

            # exon 相关：UTR5/UTR3（不做 splice 切分；与 SNVBatchRegionCount 的行为一致）
            for (es, ee) in t.exons:
                if utr5_rng[0] <= utr5_rng[1]:
                    iv = intersect_1based((es, ee), utr5_rng)
                    if iv:
                        write_iv("UTR5", chrid, iv[0], iv[1])
                if utr3_rng[0] <= utr3_rng[1]:
                    iv = intersect_1based((es, ee), utr3_rng)
                    if iv:
                        write_iv("UTR3", chrid, iv[0], iv[1])

            # CDS exon / splice / exon(body)
            cds_rng = (t.cds_start, t.cds_end)
            if cds_rng[0] <= cds_rng[1]:
                for (es, ee) in t.exons:
                    cds_iv = intersect_1based((es, ee), cds_rng)
                    if not cds_iv:
                        continue
                    cs, ce = cds_iv

                    # exonic splice edges（只对 CDS exon 生效，与 SNVBatchRegionCount 一致）
                    left = (cs, min(ce, cs + splice - 1))
                    right = (max(cs, ce - splice + 1), ce)
                    if left[0] <= left[1]:
                        write_iv("splice", chrid, left[0], left[1])
                    if right[0] <= right[1]:
                        write_iv("splice", chrid, right[0], right[1])

                    # exon body（CDS exon 去掉两端 splice）
                    body_s = cs + splice
                    body_e = ce - splice
                    if body_s <= body_e:
                        write_iv("exon", chrid, body_s, body_e)

            # intron & intronic splice flanks（对所有 exon 边界生效）
            exs = list(t.exons)
            exs.sort()
            # intronic splice flanks： [es-splice, es-1] 与 [ee+1, ee+splice]，并 clip 到 transcript body
            for (es, ee) in exs:
                a = clip_1based(es - splice, es - 1, tx_lo, tx_hi)
                b = clip_1based(ee + 1, ee + splice, tx_lo, tx_hi)
                if a:
                    write_iv("splice", chrid, a[0], a[1])
                if b:
                    write_iv("splice", chrid, b[0], b[1])

            # intron body：相邻 exon 的 gap，去掉 splice 距离
            for i in range(len(exs) - 1):
                prev_end = exs[i][1]
                next_start = exs[i+1][0]
                ins = prev_end + 1
                ine = next_start - 1
                if ins > ine:
                    continue
                body_s = ins + splice
                body_e = ine - splice
                if body_s <= body_e:
                    write_iv("intron", chrid, body_s, body_e)

    for fp in fps.values():
        fp.close()
    return paths

def run(cmd: List[str], stdin=None, stdout=None):
    p = subprocess.run(cmd, stdin=stdin, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{p.stderr}")

def bed_sort_merge(in_bed: str, out_bed: str):
    with open(out_bed, "w", encoding="utf-8") as out:
        p1 = subprocess.Popen(["bedtools", "sort", "-i", in_bed], stdout=subprocess.PIPE, text=True)
        p2 = subprocess.Popen(["bedtools", "merge", "-i", "-"], stdin=p1.stdout, stdout=out, text=True)
        p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
        if rc1 != 0 or rc2 != 0:
            raise RuntimeError("bedtools sort/merge failed")

def bed_subtract(a_bed: str, b_bed: str, out_bed: str):
    with open(out_bed, "w", encoding="utf-8") as out:
        p = subprocess.run(["bedtools", "subtract", "-a", a_bed, "-b", b_bed], stdout=out, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"bedtools subtract failed:\n{p.stderr}")

def sum_expected_overlap(gnomad_bed: str, region_bed: str) -> float:
    """
    gnomad_bed: chr start0 end0 expected
    region_bed: chr start0 end0   （已互斥化）
    返回：sum(expected * overlap_bp / window_len)
    """
    total = 0.0
    p = subprocess.Popen(
        ["bedtools", "intersect", "-a", gnomad_bed, "-b", region_bed, "-wao"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert p.stdout is not None
    for line in p.stdout:
        a = line.rstrip("\n").split("\t")
        if len(a) < 5:
            continue
        s0 = int(a[1]); e0 = int(a[2])
        exp = float(a[3])
        ov = int(a[-1])
        winlen = e0 - s0
        if winlen > 0 and ov > 0:
            total += exp * (ov / winlen)
    _, err = p.communicate()
    if p.returncode != 0:
        raise RuntimeError(f"bedtools intersect failed:\n{err}")
    return total

def read_gnomad_to_bed(in_tsv: str, out_bed: str, coord: str, expected_col: int):
    """
    生成 BED4: chr start0 end0 expected
    coord:
      - bed0: 输入 start/end 为 0-based half-open（gnomAD 常见）
      - 1based: 输入 start/end 为 1-based inclusive
    expected_col: 0-based 列号（默认示例第6列 => 5）
    """
    opener = gzip.open if in_tsv.endswith(".gz") else open
    with opener(in_tsv, "rt", encoding="utf-8", errors="ignore") as f, open(out_bed, "w", encoding="utf-8") as out:
        first = True
        for line in f:
            if not line.strip():
                continue
            if first:
                first = False
                # 跳过表头（如果第2列不是整数）
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 2 and (not is_int(parts[1])):
                    continue
            a = line.rstrip("\n").split("\t")
            if len(a) <= expected_col:
                continue
            chrom = a[0]
            if not (is_int(a[1]) and is_int(a[2])):
                continue
            start = int(a[1]); end = int(a[2])
            exp_s = a[expected_col]
            try:
                exp = float(exp_s)
            except Exception:
                continue

            # 统一 chr 前缀
            chrom = to_chr(chrom)

            if coord == "bed0":
                s0 = start
                e0 = end
            else:
                # 1-based inclusive -> bed0
                s0 = start - 1
                e0 = end
            if e0 > s0:
                out.write(f"{chrom}\t{s0}\t{e0}\t{exp}\n")

def main():
    ap = argparse.ArgumentParser(description="Summarize gnomAD window expected by gene-near regions (splice/exon/UTR/intron/up/down/intergenic) like SNVBatchRegionCount.py.")
    ap.add_argument("-i", "--infile", required=True, help="gnomAD windows TSV/TSV.GZ with columns like: chrom start end ... expected ...")
    ap.add_argument("-o", "--output", required=True, help="Output TSV: region code expected_sum")
    ap.add_argument("--gene-db", default="/data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_canonical_gene.anno.txt")
    ap.add_argument("--up", type=int, default=2000)
    ap.add_argument("--down", type=int, default=1000)
    ap.add_argument("--splice", type=int, default=2)
    ap.add_argument("--coord", choices=["bed0", "1based"], default="bed0", help="Coordinate system of infile start/end")
    ap.add_argument("--expected-col", type=int, default=5, help="0-based column index for expected (default 5 means 6th column)")
    args = ap.parse_args()

    # check bedtools
    if subprocess.run(["bash", "-lc", "command -v bedtools >/dev/null 2>&1"]).returncode != 0:
        print("ERROR: bedtools not found in PATH", file=sys.stderr)
        sys.exit(2)

    per_chr = load_transcripts(args.gene_db)

    pr_code = {
        "splice": 1,
        "exon": 2,
        "UTR5": 3,
        "UTR3": 4,
        "intron": 5,
        "upstream": 6,
        "downstream": 7,
        "intergenic": 8,
    }
    priority = ["splice","exon","UTR5","UTR3","intron","upstream","downstream"]

    with tempfile.TemporaryDirectory() as td:
        gnomad_bed = os.path.join(td, "gnomad.bed")
        read_gnomad_to_bed(args.infile, gnomad_bed, coord=args.coord, expected_col=args.expected_col)

        # raw region beds
        raw_paths = emit_raw_regions(per_chr, up=args.up, down=args.down, splice=args.splice, outdir=td)

        # sort+merge each raw bed
        merged_raw = {}
        for k, p in raw_paths.items():
            outp = os.path.join(td, f"merged_raw.{k}.bed")
            bed_sort_merge(p, outp)
            merged_raw[k] = outp

        # resolve overlaps by priority: current = raw_k - union(higher)
        resolved = {}
        higher_union = None

        for k in priority:
            cur = merged_raw[k]
            if higher_union is None:
                # first
                outp = os.path.join(td, f"resolved.{k}.bed")
                bed_sort_merge(cur, outp)
                resolved[k] = outp
                higher_union = outp
            else:
                tmp_sub = os.path.join(td, f"tmp_sub.{k}.bed")
                bed_subtract(cur, higher_union, tmp_sub)
                outp = os.path.join(td, f"resolved.{k}.bed")
                bed_sort_merge(tmp_sub, outp)
                resolved[k] = outp

                # update higher_union = merge(higher_union + outp)
                catp = os.path.join(td, f"tmp_cat_union.{k}.bed")
                with open(catp, "w", encoding="utf-8") as out:
                    with open(higher_union, "r", encoding="utf-8") as f1:
                        out.write(f1.read())
                    with open(outp, "r", encoding="utf-8") as f2:
                        out.write(f2.read())
                new_union = os.path.join(td, f"union_upto.{k}.bed")
                bed_sort_merge(catp, new_union)
                higher_union = new_union

        # sum expected by region (disjoint by construction)
        sums = {}
        total_expected = 0.0
        with open(gnomad_bed, "r", encoding="utf-8") as f:
            for line in f:
                a = line.rstrip("\n").split("\t")
                if len(a) >= 4:
                    total_expected += float(a[3])

        used = 0.0
        for k in priority:
            s = sum_expected_overlap(gnomad_bed, resolved[k])
            sums[k] = s
            used += s

        intergenic = max(0.0, total_expected - used)
        sums["intergenic"] = intergenic

        # write output in code order
        order = ["splice","exon","UTR5","UTR3","intron","upstream","downstream","intergenic"]
        with open(args.output, "w", encoding="utf-8") as out:
            out.write("region\tcode\texpected_sum\n")
            for k in order:
                out.write(f"{k}\t{pr_code[k]}\t{(sums.get(k, 0.0)):.10f}\n")

if __name__ == "__main__":
    main()