import sys

########Functions############################################################
# transformation for minus strand
def toTrans(x):
    x=x.upper()
    x=list(x)
    for i in range(len(x)):
        if x[i]=='A':
            x[i]='T'
        elif x[i]=='T':
            x[i]='A'        
        elif x[i]=='C':
            x[i]='G'
        elif x[i]=='G':
            x[i]='C'
    x.reverse()
    return(''.join(x))
# 根据SNV突变位点获取改变后的序列
def getSNVSeq(seq,pos_real,pos,alt):
    seq=list(seq)
    pos=int(pos)
    seq[pos_real.index(pos)]=alt
    return(''.join(seq))
# 根据MNV突变位点获取改变后的序列
def getMNVSeq(seq,pos_real,pos_list,alt_list):
    seq=list(seq)
    pos_list=[int(i) for i in pos_list]
    for i in range(len(pos_list)):
        seq[pos_real.index(pos_list[i])]=alt_list[i]
    return(''.join(seq))
####################################################################
# python 002getSeq_total_UTR3.py ../01data/en_GRCh38.108_gene.anno.txt ../01data/chr_fix.fa ./ &

region_url = sys.argv[1] # region_url='../01data/en_GRCh38.108_gene.anno.txt'
fa_url = sys.argv[2] # fa_url='../01data/chr_fix.fa'
output_file = sys.argv[3] # output_file='./'


####################################################################

########Data loading, processing and output############################################################
# Load Sequence
f=open(fa_url)
s=f.readlines()
chr_seq = {}
for i in range(0,len(s),2):
    chr_seq[s[i].strip('\n').replace('>','')]=s[i+1].strip('\n')
del s
f.close()


# Obtain sequences
with open(region_url) as f, open(output_file+'tar_ref.seq','w') as result1, open(output_file+'miranda_ref.seq','w') as result4:
    # Transcript ID (unique),Chromosome name ,Strand ,Gene ID ,Gene symbol ,Start (order) ,End ,CDS Start ,CDS End ,Number of exons ,Starts of exons ,Ends of exons ,(Start-2000)-(End+1000) ,Upstream ,Downstream ,5'UTR ,3'UTR ,Exons ,Splices
    # ENST00000361739 MT + ENSG00000198712 MT-CO2  7586    8269    7586 8269 1 7586 8269 5586-9269 5586-7585 8270-9269 1-7586-8269
    for i in f:
        line=i.strip('\n').split('\t')
        chrid=line[1]
        ncRNA_ID=line[0]
        ncRNA_strand=line[2]
        # 如果没有UTR3，跳过
        if line[16]=='':
            continue
        exons=line[16].split(',')
        # 提取序列
        if ncRNA_strand=='+':
            # exons=['1-10-20','2-40-60']
            seq=''
            pos_real=[]
            for ii in exons:
                seq=seq+chr_seq[chrid][(int(ii.split('-')[1])-1):int(ii.split('-')[2])]
            for ii in exons:
                pos_real=pos_real+list(range(int(ii.split('-')[1]),int(ii.split('-')[2])+1))
        else:
            # exons=['1-20-15','2-10-4'] -
            seq=[]
            pos_real=[]
            for ii in exons:
                seq.append(chr_seq[chrid][(int(ii.split('-')[2])-1):int(ii.split('-')[1])])
            seq.reverse()
            seq=''.join(seq)
            for ii in exons:
                pos_real=pos_real+list(range(int(ii.split('-')[2]),int(ii.split('-')[1])+1))
            pos_real.sort()        

        seq=seq.upper()
        # 获取ref的seq
        ref=seq
        if ncRNA_strand == '-':
            ref = toTrans(ref)

        result1.write('|'.join([ncRNA_ID,ncRNA_strand])+'\t9606\t'+ref.replace("T", "U")+'\n')
        result4.write('>'+'|'.join([ncRNA_ID,ncRNA_strand])+'\n')
        result4.write(ref.replace("T", "U")+'\n')
####################################################################
