## Mutation mechanism selection whole genome sequencing data (1000G and GTEx)

## Data preprocessing
## Extract vcf snv data and AC data # Get 1000G_snv.txt GTEx_snv.txt
    python 01getSNVAC.py /data/jinww/04reference/publicDB/ /data/jinww/mnv/analyse3/03basePattern/01process/
## Split from 1000G.txt and GTEx.txt to 1000G/ and GTEx/, and match SNV AC data to MNV
    python 02splitData.py 1000G.txt &
    python 02splitData.py GTEx.txt &
## Merge the two datasets
    python 03merge.py 1000G_adjust.txt GTEx_adjust.txt mnv.txt
## Extract sequence ±5bp append (10bp)
    python 04addSeq.py /data/jinww/mnv/library/human/01gtf/ensembl/hg38/chr_fix.fa 5 mnv.txt mnv_5.txt &
