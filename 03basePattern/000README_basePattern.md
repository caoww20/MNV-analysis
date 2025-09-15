## Mutation mechanism selection whole genome sequencing data (1000G and GTEx)

## Data preprocessing
## Link 1000G.txt and GTEx.txt files to part3
ln -s /data/jinww/mnv/analyse3/02data/02datasets/GTEx.txt ./
ln -s /data/jinww/mnv/analyse3/02data/02datasets/1000G.txt ./
## Extract vcf snv data and AC data # Get 1000G_snv.txt GTEx_snv.txt
    python 01getSNVAC.py /data/jinww/04reference/publicDB/ /data/jinww/mnv/analyse3/03basePattern/01process/
## Split from 1000G.txt and GTEx.txt to 1000G/ and GTEx/, and match SNV AC data to MNV (here filter out multi-allelic, because theoretically multi-allelic cannot be one-step, but theoretically possible)
    python 02splitData.py 1000G.txt &
    python 02splitData.py GTEx.txt &
## Merge the two datasets
    python 03merge.py 1000G_adjust.txt GTEx_adjust.txt mnv.txt
## Extract sequence ±5bp append (10bp)
    python 04addSeq.py /data/jinww/mnv/library/human/01gtf/ensembl/hg38/chr_fix.fa 5 mnv.txt mnv_5.txt &
## Judge one-step and repeats
    python 05judgeMec.py mnv_5.txt &
    awk '$14=="1" {print $o}' mnv_5_adjust.txt |wc -l # Check if the result meets expectations


## Global pattern statistics
    Conclusion: 2-joint MNVs are significantly higher at distance 1 than other distances, others have delayed peaks but not significant, partly due to distance limitation in identification, partly because there may not be such a prominent mechanism as 2-joint
## Statistics of adjacent 2-joint patterns (including mutation patterns per unit)
    Adjacent MNV trends are similar to reports, for non-adjacent MNVs, ANA->GNG, CNC->TNT and their reverse complements account for the highest, followed by transitions.
## Correlation of curves for different joints after excluding adjacent 2-joints, subsequent analysis does not target adjacent 2-joints
## Further analysis of one-step repeat junction (per joint unit, per distance unit): includes both cross-sectional and longitudinal
## Further analysis of ATCG percentages (per joint unit)


## Calculate LOF, using Luo Haohui's data
cp /home/luohh/Jinww/Loffee 02LOF
