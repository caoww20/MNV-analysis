# base:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/merge/gwas.txt
# gene:/data/jinww/mnv/analyse3/02data/02annotation/gene_up.txt
# gwas:/data/jinww/04reference/publicDB/gwas/download/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv
# NCF:/data/jinww/04reference/publicDB/cancer_gene/NCG_cancerdrivers_annotation_supporting_evidence.tsv
# intOGen:/data/jinww/04reference/publicDB/cancer_gene/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv
# 012:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/012/11Cancers.survival.p.txt
# 01:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/01/11Cancers.survival.p.txt
# 02:/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/02/11Cancers.survival.p.txt
# python merge.py

# 读入base数据中的MNVid
mnvid={}
mnvidCancer={}
with open('/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/merge/gwas.txt','r') as f:
    # 跳过第一行
    next(f)
    for line in f:
        a=line.strip().split('\t')
        mnvid[a[1]]=[]
        mnvidCancer[a[0]+'|'+a[1]]=['NA','NA','NA']

# 读入gene数据中的MNVid
# with open('/data/jinww/mnv/analyse3/02data/02annotation/gene_up.txt','r') as f: # hg38
with open('/data/jinww/mnv/database/hg19/anno/gene.txt','r') as f:
    for line in f:
        a=line.strip().split('\t')
        if a[2] in mnvid:
            rsid=a[5]
            if a[-1] != '.':
                b=a[-1].split('|')
                anno=[]
                for i in b:
                    c=i.split(' ')
                    gene_symbol=c[1]
                    location=c[3]
                    anno.append(gene_symbol+'-'+location)
                # anno去重
                anno=list(set(anno))
                anno='&'.join(anno)
                mnvid[a[2]]=[rsid+'|'+anno]
            else:
                mnvid[a[2]]=[rsid+'|'+'intergenic']
# 从mnvid中获取rsid的集合
rsid={}
for i in mnvid:
    a=mnvid[i][0].split('|')[0].split(',')
    for j in a:
        if j !='.':
            rsid[j]=0
# 确定rsid是否存在于gwas数据中
with open('/data/jinww/04reference/publicDB/gwas/download/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv','r') as f:
    for line in f:
        # n=line.strip().split('\t')
        a=line.strip().split('\t')[21]
        if a in rsid:
            rsid[a]=1
# 从mnvid中获取基因的集合
gene={}
for i in mnvid:
    a=mnvid[i][0].split('|')[1]
    if a !='intergenic':
        b=a.split('&')
        for j in b:
            c=j.split('-')
            gene[c[0]]=[0,0]
# 根据NCF的数据，确定gene是否是cancer gene
with open('/data/jinww/04reference/publicDB/cancer_gene/NCG_cancerdrivers_annotation_supporting_evidence.tsv','r') as f:
    next(f)
    for line in f:
        a=line.strip().split('\t')
        if a[1] in gene:
            gene[a[1]][0]=1
# 根据intOGen的数据，确定gene是否是cancer gene
with open('/data/jinww/04reference/publicDB/cancer_gene/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv','r') as f:
    next(f)
    for line in f:
        a=line.strip().split('\t')
        if a[0] in gene:
            gene[a[0]][1]=1
# 读入012survival数据
with open('/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/012/11Cancers.survival.p.txt','r') as f:
    next(f)
    for line in f:
        a=line.strip().split(' ')
        id=a[0]+'|'+a[1]
        if id in mnvidCancer:
            mnvidCancer[id][0]=a[2]
# 读入01survival数据
with open('/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/01/11Cancers.survival.p.txt','r') as f:
    next(f)
    for line in f:
        a=line.strip().split(' ')
        id=a[0]+'|'+a[1]
        if id in mnvidCancer:
            mnvidCancer[id][1]=a[2]
# 读入02survival数据
with open('/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/MNV/04.Plot/01.SurvivalPlot/10years/02/11Cancers.survival.p.txt','r') as f:
    next(f)
    for line in f:
        a=line.strip().split(' ')
        id=a[0]+'|'+a[1]
        if id in mnvidCancer:
            mnvidCancer[id][2]=a[2]
# 读入base数据
with open('/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/merge/gwas.txt','r') as f,open('/data/jinww/mnv/analyse3/08cancer/UKB50w/01newdata/merge/res_hg19.txt','w') as out:
    out.write('CancerType\tMNVID\tRSID\tAF\tBeta\tP(GWAS)\tP012(KM)\tP01(KM)\tP02(KM)\tGene\tLocation\tCancerGene(NCG)\tCancerGene(intOGen)\tinGWASCatelog\n')
    next(f)
    for line in f:
        a=line.strip().split('\t')
        mnvid2=a[1]
        id=a[0]+'|'+a[1]
        rsid2=mnvid[mnvid2][0].split('|')[0]
        gene_location=mnvid[mnvid2][0].split('|')[1]
        survival=mnvidCancer[id]
        inGWASCatelog=0
        for i in rsid2.split(','):
            if i !='.':
                if rsid[i] == 1:
                    inGWASCatelog=1
                    break
        if gene_location !='intergenic':
            for i in gene_location.split('&'):
                # NEK10-up:518,up:526&NEK10-intron,intron
                iscancerGene=[0,0]
                gene_symbol=i.split('-')[0]
                location=i.split('-')[1]
                if gene_symbol in gene:
                    iscancerGene= gene[gene_symbol]
                newa=[a[0],mnvid2,rsid2,a[4],a[2],a[3],survival[0],survival[1],survival[2],gene_symbol,location,iscancerGene[0],iscancerGene[1],inGWASCatelog]
                # 将newa转化为str
                newa=[str(i) for i in newa]
                out.write('\t'.join(newa)+'\n')
        else:
            newa=[a[0],mnvid2,rsid2,a[4],a[2],a[3],survival[0],survival[1],survival[2],'-','-','-','-',inGWASCatelog]
            newa=[str(i) for i in newa]
            out.write('\t'.join(newa)+'\n')
                
        