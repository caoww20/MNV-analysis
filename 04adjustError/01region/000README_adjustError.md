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

## Calculate conservation of MNV sites (this part does not need to be updated, using old data)
## Run getSNPfromHG38.py to obtain positions of each SNP
    # total
    python 001getTotalSNP.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv ./02pos/total.pos &
    # gene 
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/gene.txt gene exon ./02pos/cds.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/gene.txt gene intro ./02pos/intro.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/gene.txt gene splice ./02pos/splice.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/gene.txt gene UTR3 ./02pos/utr3.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/gene.txt gene UTR5 ./02pos/utr5.pos &
    # non 
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/circ.txt other 0 ./02pos/circ.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/lnc.txt other 0 ./02pos/lnc.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/pre_miRNA.txt other 0 ./02pos/miRNA.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/piR.txt other 0 ./02pos/piR.pos &
    # reg
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/enhancer.txt other 0 ./02pos/enhancer.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/atac.txt other 0 ./02pos/atac.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/CE.txt other 0 ./02pos/CE.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/miRBS.txt other 0 ./02pos/miRBS.pos &
    python 002getSNPfromHG38.py /data/jinww/mnv/analyse3/02data/02annotation/TFBS.txt other 0 ./02pos/TFBS.pos &
# Use linux sort to deduplicate each
    ls *pos|while read id;do (sort -u $id >$id.sort &);done
# Get intergenic
    cat cds.pos.sort splice.pos.sort utr3.pos.sort utr5.pos.sort intro.pos.sort lnc.pos.sort miRNA.pos.sort circ.pos.sort piR.pos.sort  |sort -u >gene.txt
    python 003getIntergenic.py ./02pos/total.pos.sort ./02pos/gene.txt ./02pos/intergenic.pos
    sort -u intergenic.pos >intergenic.pos.sort &
## Run 04getCos.py
    for i in {1..22}; do ( python 004getCos.py 01SNPConservationFromUCSC/chr${i}.wig $i 02pos/total.pos.sort 03pos_cos/total_${i} &); done
    cat total* >cos.score
    rm total*
## Run 05getRegionCos.py to get data for each region
    ls  *sort|while read id;do ( python ../005getRegionCos.py cos.score $id ${id%.*}.cos &);done
    cat *cos >res.cos
## Results show that circRNA results are somewhat abnormal, this is related to the length of circ itself (average of circRNAs within 10KB is 0.097, but decreases gradually with length increase), generally can reflect conservation, the stronger the conservation
    for i in {1..11}; do (bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i /data/jinww/mnv/analyse3/02data/01origin/hg38mnv -t non -D hg38_circ_aj_${i}.txt -o circ_${i}.txt &);done &
    for i in {1..11}; do (python 002getSNPfromHG38.py 04circ_diff/circ_${i}.txt other circRNA 04circ_diff/circ_${i}.pos &);done &
    ls *pos|while read id;do (sort -u $id >$id.sort &);done
    ls  04circ_diff/*sort|while read id;do ( python 005getRegionCos.py 03pos_cos/cos.score $id ${id%.*}.cos &);done


