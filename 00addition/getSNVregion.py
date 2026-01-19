import sys
import re

# 新用法：
#   python 002getSummary.py <snv_file> <annotation_dir/> <out_res>
#
# snv_file: 每行至少包含一个 SNV 主键，如 "13:16000187:C:T"
# annotation_dir: 包含 gene.txt / lnc.txt / pre_miRNA.txt / circ.txt / piR.txt / atac.txt / CE.txt / enhancer.txt / miRBS.txt / TFBS.txt
# out_res: 输出统计结果

snv_file = sys.argv[1]
url = sys.argv[2]
res = sys.argv[3]

if not url.endswith("/"):
    url += "/"

# 匹配常见 SNV key：chr:pos:ref:alt（chr允许 1-22,X,Y,MT 或带 chr 前缀）
SNV_TOKEN_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT):([0-9]+):([ACGTN]+):([ACGTN]+)$", re.IGNORECASE)

def normalize_snv(token: str) -> str | None:
    """把各种写法统一成 'chr:pos:ref:alt'（不带 chr 前缀），返回 None 表示不是 SNV 格式。"""
    token = token.strip()
    m = SNV_TOKEN_RE.match(token)
    if not m:
        return None
    chrom, pos, ref, alt = m.group(1).upper(), m.group(2), m.group(3).upper(), m.group(4).upper()
    # 统一不加 'chr'：比如 13:16000187:C:T
    return f"{chrom}:{pos}:{ref}:{alt}"

def extract_snv_from_fields(fields: list[str]) -> str | None:
    """从一行拆分字段中提取 SNV key。优先寻找符合 chr:pos:ref:alt 的字段。"""
    for x in fields:
        k = normalize_snv(x)
        if k:
            return k
    return None

def load_snv_set(path: str) -> set[str]:
    snvs = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            k = extract_snv_from_fields(fields)
            if k:
                snvs.add(k)
    return snvs

def create_anno():
    return {
        "UTR5": 0,
        "exon": 0,
        "splice": 0,
        "intron": 0,
        "UTR3": 0,
        "lncRNA": 0,
        "miRNA": 0,
        "circRNA": 0,
        "piRNA": 0,
        "ATAC": 0,
        "conserved_element": 0,
        "enhancer": 0,
        "miRBS": 0,
        "TFBS": 0,
        "intergenic": 0,
    }

snv_set = load_snv_set(snv_file)

# 收集各区域 SNV（用 set 去重）
UTR5, exon, splice, intron, UTR3 = set(), set(), set(), set(), set()
lnc, miR, circRNA, piRNA = set(), set(), set(), set()
ATAC, conserved_element, enhancer, miRBS, TFBS = set(), set(), set(), set(), set()

def add_if_hit(path: str, predicate, bucket: set[str]):
    """通用读取：若该行能提取到 snv_key 且在 snv_set 中，并且 predicate(fields) 为真，则加入 bucket。"""
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # 允许注释文件是空格分隔时退化
            if len(fields) == 1:
                fields = line.split()

            snv_key = extract_snv_from_fields(fields)

            # 兼容旧注释格式：如果找不到 token，就尝试用第3列（a[2]）当作 key
            if snv_key is None and len(fields) > 2:
                snv_key = normalize_snv(fields[2])  # 如果 a[2] 恰好是 13:pos:ref:alt
            if snv_key is None:
                continue

            if snv_key in snv_set and predicate(fields):
                bucket.add(snv_key)

# gene.txt：按 consequence 分类
gene_path = url + "gene.txt"
with open(gene_path) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) == 1:
            fields = line.split()

        snv_key = extract_snv_from_fields(fields)
        if snv_key is None and len(fields) > 2:
            snv_key = normalize_snv(fields[2])
        if snv_key is None or snv_key not in snv_set:
            continue

        cons = fields[-1] if fields else ""
        if cons == ".":
            continue

        if "UTR5_variant" in cons:
            UTR5.add(snv_key)
        if "exon" in cons:
            exon.add(snv_key)
        if "splice_" in cons:
            splice.add(snv_key)
        if "intron_variant" in cons:
            intron.add(snv_key)
        if "UTR3_variant" in cons:
            UTR3.add(snv_key)

# 其他功能文件：按最后一列标签判断
add_if_hit(url + "lnc.txt", lambda a: ("all lncRNA" in a[-1]), lnc)
add_if_hit(url + "pre_miRNA.txt", lambda a: ("all miRNA" in a[-1]), miR)
add_if_hit(url + "circ.txt", lambda a: ("all circRNA" in a[-1]), circRNA)
add_if_hit(url + "piR.txt", lambda a: ("all piRNA" in a[-1]), piRNA)
add_if_hit(url + "atac.txt", lambda a: ("all ATAC" in a[-1]), ATAC)
add_if_hit(url + "CE.txt", lambda a: ("all conserved_element" in a[-1]), conserved_element)
add_if_hit(url + "enhancer.txt", lambda a: ("all enhancer" in a[-1]), enhancer)
add_if_hit(url + "miRBS.txt", lambda a: ("all miRBS" in a[-1]), miRBS)
add_if_hit(url + "TFBS.txt", lambda a: ("all TFBS" in a[-1]), TFBS)

anno = create_anno()
anno["UTR5"] = len(UTR5)
anno["exon"] = len(exon)
anno["splice"] = len(splice)
anno["intron"] = len(intron)
anno["UTR3"] = len(UTR3)
anno["lncRNA"] = len(lnc)
anno["miRNA"] = len(miR)
anno["circRNA"] = len(circRNA)
anno["piRNA"] = len(piRNA)
anno["ATAC"] = len(ATAC)
anno["conserved_element"] = len(conserved_element)
anno["enhancer"] = len(enhancer)
anno["miRBS"] = len(miRBS)
anno["TFBS"] = len(TFBS)

# intergenic：不在任何已标注集合里的 SNV
annotated = set().union(
    UTR5, exon, splice, intron, UTR3,
    lnc, miR, circRNA, piRNA,
    ATAC, conserved_element, enhancer, miRBS, TFBS
)
intergenic = snv_set - annotated
anno["intergenic"] = len(intergenic)

with open(res, "w") as out:
    for k, v in anno.items():
        out.write(f"{k}\t{v}\n")