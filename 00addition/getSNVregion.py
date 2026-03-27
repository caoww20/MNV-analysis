import sys
import re

# New usage:
#   python 002getSummary.py <snv_file> <annotation_dir/> <out_res>
#
# snv_file: each line contains at least one SNV key, e.g. "13:16000187:C:T"
# annotation_dir: contains gene.txt / lnc.txt / pre_miRNA.txt / circ.txt / piR.txt / atac.txt / CE.txt / enhancer.txt / miRBS.txt / TFBS.txt
# out_res: output stats

snv_file = sys.argv[1]
url = sys.argv[2]
res = sys.argv[3]

if not url.endswith("/"):
    url += "/"

# Match common SNV key: chr:pos:ref:alt (chr allows 1-22,X,Y,MT or with chr prefix)
SNV_TOKEN_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT):([0-9]+):([ACGTN]+):([ACGTN]+)$", re.IGNORECASE)

def normalize_snv(token: str) -> str | None:
    """Normalize variants to 'chr:pos:ref:alt' (no chr prefix); return None if not SNV format."""
    token = token.strip()
    m = SNV_TOKEN_RE.match(token)
    if not m:
        return None
    chrom, pos, ref, alt = m.group(1).upper(), m.group(2), m.group(3).upper(), m.group(4).upper()
    # Normalize without 'chr', e.g. 13:16000187:C:T
    return f"{chrom}:{pos}:{ref}:{alt}"

def extract_snv_from_fields(fields: list[str]) -> str | None:
    """Extract SNV key from split fields, preferring chr:pos:ref:alt tokens."""
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

# Collect SNVs per region (use set for de-duplication)
UTR5, exon, splice, intron, UTR3 = set(), set(), set(), set(), set()
lnc, miR, circRNA, piRNA = set(), set(), set(), set()
ATAC, conserved_element, enhancer, miRBS, TFBS = set(), set(), set(), set(), set()

def add_if_hit(path: str, predicate, bucket: set[str]):
    """Generic loader: add snv_key to bucket if in snv_set and predicate(fields) is True."""
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # Allow fallback if annotation file is space-delimited
            if len(fields) == 1:
                fields = line.split()

            snv_key = extract_snv_from_fields(fields)

            # Backward-compatible annotation format: fall back to column 3 as key
            if snv_key is None and len(fields) > 2:
                snv_key = normalize_snv(fields[2])  # If fields[2] matches 13:pos:ref:alt
            if snv_key is None:
                continue

            if snv_key in snv_set and predicate(fields):
                bucket.add(snv_key)

# gene.txt: classify by consequence
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

# Other annotation files: check the last column tag
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

# intergenic: SNVs not in any annotated set
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