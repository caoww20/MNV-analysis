import sys
file=sys.argv[1] # '/data/jinww/mnv/analyse3/02data/01origin/hg38mnv'
url=sys.argv[2] # '/data/jinww/mnv/analyse3/02data/02annotation/'
dataset=sys.argv[3]
maf_switch=int(sys.argv[4]) # 如果是1的话，启用MAF过滤
res=sys.argv[5]
dataset=dataset.split(',')
        
# 读入humanMNV的数据库
mnv=[]
with open(file) as f:
    for i in f:
        a=i.strip().split('\t')
        b=a[-1]
        # if 'TCGA' in b or 'GTEx' in b or '1000G' in b or 'UKB20w' in b or 'UKB50w' in b:
        #   mnv.append(a[2])
        flag=0
        for ii in dataset:
            if ii in b:
                flag=1
                break
        
        if maf_switch==1:
            if flag==1:
                for ii in b.split('|'):
                    if ii.startswith(dataset[0]):
                        ac=int(ii.split('/')[0].split(':')[1])
                        break
                if ac>3: # 相当于maf=0.05%
                    mnv.append(a[2])
        else:
            if flag==1:
                mnv.append(a[2])

mnv=set(mnv)

# protein
def createGeneAnno():
    anno={}
    anno['UTR5'] = 0
    anno['exon'] = 0
    anno['splice'] = 0
    anno['intron'] = 0
    anno['UTR3'] = 0
    return anno
anno=createGeneAnno()

UTR5 = []
exon = []
splice = []
intron = []
UTR3 = []

with open(url+'gene.txt') as f:
    for i in f:
        a=i.strip().split('\t')
        if a[-1] !='.' and a[2] in mnv:
            if 'UTR5_variant' in a[-1]:
                UTR5.append(a[2])
            if 'exon' in a[-1]:
                exon.append(a[2])
            if 'splice_' in a[-1]:
                splice.append(a[2])
            if 'intron_variant' in a[-1]:
                intron.append(a[2])
            if 'UTR3_variant' in a[-1]:
                UTR3.append(a[2])                                     


UTR5 = set(UTR5)
exon = set(exon)
splice = set(splice)
intron = set(intron)
UTR3 = set(UTR3)

anno['UTR5']=len(UTR5)
anno['exon']=len(exon)
anno['splice']=len(splice)
anno['intron']=len(intron)
anno['UTR3']=len(UTR3)

# lncRNA
lnc=[]
f=open(url+'lnc.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all lncRNA' in a[-1]:
        lnc.append(a[2])
f.close()
lnc=set(lnc)
anno['lncRNA']=len(lnc)

# miR 
miR=[]
f=open(url+'pre_miRNA.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all miRNA' in a[-1]:
        miR.append(a[2])
f.close()
miR=set(miR)
anno['miRNA']=len(miR)

# circRNA
circRNA=[]
f=open(url+'circ.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all circRNA' in a[-1]:
        circRNA.append(a[2])
f.close()
circRNA=set(circRNA)
anno['circRNA']=len(circRNA)

# piRNA
piRNA=[]
f=open(url+'piR.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all piRNA' in a[-1]:
        piRNA.append(a[2])
f.close()
piRNA=set(piRNA)
anno['piRNA']=len(piRNA)

# ATAC
ATAC=[]
f=open(url+'atac.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all ATAC' in a[-1]:
        ATAC.append(a[2])
f.close()
ATAC=set(ATAC)
anno['ATAC']=len(ATAC)

# CE
conserved_element=[]
f=open(url+'CE.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all conserved_element' in a[-1]:
        conserved_element.append(a[2])
f.close()
conserved_element=set(conserved_element)
anno['conserved_element']=len(conserved_element)

# enhancer
enhancer=[]
f=open(url+'enhancer.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all enhancer' in a[-1]:
        enhancer.append(a[2])
f.close()
enhancer=set(enhancer)
anno['enhancer']=len(enhancer)

# miRBS
miRBS=[]
f=open(url+'miRBS.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all miRBS' in a[-1]:
        miRBS.append(a[2])
f.close()
miRBS=set(miRBS)
anno['miRBS']=len(miRBS)

# TFBS
TFBS=[]
f=open(url+'TFBS.txt')
for i in f:
    a=i.strip().split('\t')
    if a[2] in mnv and 'all TFBS' in a[-1]:
        TFBS.append(a[2])
f.close()
TFBS=set(TFBS)
anno['TFBS']=len(TFBS)

# intergenic
# total - protein -lncRNA -miRNA -circRNA - piRNA 
intergenic=mnv-UTR5-exon-splice-intron-UTR3-lnc-miR-circRNA-piRNA
anno['intergenic']=len(intergenic)


f=open(res,'w')
for i in anno:
    f.write(i+'\t'+str(anno[i])+'\n')
f.close()