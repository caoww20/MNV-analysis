import sys
file_url=sys.argv[1] # '/data/jinww/mnv/MNVAnno/database/human/'
res=sys.argv[2] # region.txt
####################################################################
## 获取每个区域的基因组覆盖度从而计算MNV的覆盖度 MNV数量/基因组区域
def merge_intervals(intervals):
    if len(intervals) == 0:
        return []
    intervals.sort(key=lambda x: x[0])
    result = [intervals[0]]
    for i in range(1, len(intervals)):
        if result[-1][1] >= intervals[i][0]:
            result[-1][1] = max(result[-1][1], intervals[i][1])
        else:
            result.append(intervals[i])
    return result
# intervals = [[5,10],[15,18],[1,3],[2,6]]
# print(merge_intervals(intervals))
def sum_range(lst):
    diff_sum = 0
    for sublist in lst:
        a, b = sublist
        diff_sum += (b - a + 1) # 因为是左闭右闭
    return(diff_sum)
# lst = [[1, 5], [10, 20]]
# sum_range(lst)
##intergenic
def getintergenic(url):
    chr_list={}
    for i in range(1,23,1):
        chr_list[str(i)]=[]
    # gene
    f=open(url+'en_GRCh38.108_gene.anno.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            chr_list[chrid].append([int(a[5]),int(a[6])])
    f.close()
    # lnc
    f=open(url+'en_GRCh38.108_non.anno.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]=='lncRNA':
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()
    f=open(url+'hg38_lnc_aj.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]=='lncRNA':
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()
    # miR
    f=open(url+'en_GRCh38.108_non.anno.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]=='miRNA':
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()
    f=open(url+'hg38_miRNA_aj.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]=='miRNA':
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()
    # circ
    f=open(url+'hg38_circ_aj.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]=='circRNA':
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()
    # piRNA
    f=open(url+'hg38_piR_aj.txt')
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]=='piRNA':
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()

    for i in chr_list:
        chr_list[i]=merge_intervals(chr_list[i])
    for i in chr_list:
        chr_list[i]=sum_range(chr_list[i])
    
    all_sum=0
    for i in chr_list:
        all_sum +=  chr_list[i]
    return(all_sum)
# other region
def getbed(input1,flag,dataNum,input2):
    chr_list={}
    for i in range(1,23,1):
        chr_list[str(i)]=[]
    # chr_list['X']=[]
    f=open(input1)
    for i in f:
        a=i.strip().split('\t')
        chrid=a[1]
        if chrid in chr_list:
            if a[0]==flag:
                chr_list[chrid].append([int(a[2]),int(a[3])])
    f.close()
    if dataNum == 2:
        f=open(input2)
        for i in f:
            a=i.strip().split('\t')
            chrid=a[1]
            if chrid in chr_list:
                if a[0]==flag:
                    chr_list[chrid].append([int(a[2]),int(a[3])])
        f.close()
    for i in chr_list:
        chr_list[i]=merge_intervals(chr_list[i])
    for i in chr_list:
        chr_list[i]=sum_range(chr_list[i])
    all_sum=0
    for i in chr_list:
        all_sum +=  chr_list[i]
    return(all_sum)
##gene

def getPart(k):
    chr_list={}
    for i in range(1,23,1):
        chr_list[str(i)]=[]
    f=open(file_url+'en_GRCh38.108_gene.anno.txt')
    for i in f:
        a=i.strip().split('\t')
        if a[k] !='' and a[1] in chr_list:
            for ii in a[k].split(','):
                b=ii.split('-')
                if len(b)==3:
                    tmp=[int(b[1]),int(b[2])]
                    tmp.sort()
                    chr_list[a[1]].append(tmp)
                else:
                    tmp=[int(b[1]),int(b[1])]
                    tmp.sort()
                    chr_list[a[1]].append(tmp)
                    tmp=[int(b[3]),int(b[4])]
                    tmp.sort()
                    chr_list[a[1]].append(tmp)
    for i in chr_list:
        chr_list[i]=merge_intervals(chr_list[i])
    for i in chr_list:
        chr_list[i]=sum_range(chr_list[i])
    all_sum=0
    for i in chr_list:
        all_sum +=  chr_list[i]
    return(all_sum)
utr5=getPart(15)
utr3=getPart(16)
cds=getPart(17)
splice=getPart(18)

def getTrans():
    chr_list={}
    for i in range(1,23,1):
        chr_list[str(i)]=[]
    f=open(file_url+'en_GRCh38.108_gene.anno.txt')
    for i in f:
        a=i.strip().split('\t')
        if a[1] in chr_list:
            chr_list[a[1]].append([int(a[5]),int(a[6])])
    f.close()
    for i in chr_list:
        chr_list[i]=merge_intervals(chr_list[i])
    for i in chr_list:
        chr_list[i]=sum_range(chr_list[i])
    trans=0
    for i in chr_list:
        trans +=  chr_list[i]
    return(trans)
trans=getTrans()

intro=trans-utr5-cds-utr3-splice


##lncRNA
lncRNA=getbed(file_url+'en_GRCh38.108_non.anno.txt','lncRNA',2,file_url+'hg38_lnc_aj.txt')
##miRNA
miRNA=getbed(file_url+'en_GRCh38.108_non.anno.txt','miRNA',2,file_url+'hg38_miRNA_aj.txt')
##circRNA
circRNA=getbed(file_url+'hg38_circ_aj.txt','circRNA',1,'')
##piRNA
piRNA=getbed(file_url+'hg38_piR_aj.txt','piRNA',1,'')


##ATAC
ATAC=getbed(file_url+'hg38_ATAC_aj.txt','ATAC',1,'')
##CE
CE=getbed(file_url+'hg38_CE_aj.txt','conserved_element',1,'')
##enhancer
enhancer=getbed(file_url+'hg38_enhancer_aj.txt','enhancer',1,'')
##miRBS
miRBS=getbed(file_url+'hg38_miRBS_aj.txt','miRBS',1,'')
##TFBS
TFBS=getbed(file_url+'hg38_TFBS_aj.txt','TFBS',1,'')


## intergenic
intergenic=2875001522-getintergenic(file_url)
# intergenic=0

regionNames=['UTR5','cds','splice','intron','UTR3','lncRNA','miRNA','circRNA','piRNA','ATAC','conserved_element','enhancer','miRBS','TFBS','intergenic']
regionBases=[utr5,cds,splice,intro,utr3,lncRNA,miRNA,circRNA,piRNA,ATAC,CE,enhancer,miRBS,TFBS,intergenic]


# 读入数量
f=open(res,'w')
for i in range(len(regionBases)):
    f.write(regionNames[i]+'\t'+str(regionBases[i])+'\n')
f.close()




