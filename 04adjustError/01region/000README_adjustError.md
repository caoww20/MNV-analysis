## Genome region distribution =========================
## Run 01getRegion.py to obtain lengths of different regions
    python 001getRegion.py /data/jinww/mnv/MNVAnno/database/human/ region.txt
## Run 02getSummary.py to get MNV counts (here only extract data that completely falls within regions, all tags)
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ 1000G,GTEx,TCGA,UKB20w,UKB50w 0 ./01num/all.txt &
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ GnomAD_WGS,GnomAD_WES 0 ./01num/GnomAD.txt &
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ GnomAD_WGS 0 ./01num/gnomAD_WGS.txt &
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ GTEx 0 ./01num/GTEx.txt &
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ TCGA 0 ./01num/TCGA.txt &
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ UKB20w 0 ./01num/UKB20w.txt &
    python 002getSummary.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv /data/jinww/mnv/analyse3/02data/02annotation/ UKB50w 0 ./01num/UKB50w.txt &
