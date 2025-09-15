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
# python 002getSeq_total_lncRNA.py /data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt ../01data/chr_fix.fa ./ &

region_url = sys.argv[1] # region_url='/data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt'
fa_url = sys.argv[2] # fa_url='chr_fix.fa'
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


# Obtain 5'utr and cds sequences
with open(region_url) as f, open(output_file+'tar_ref.seq','w') as result1, open(output_file+'miranda_ref.seq','w') as result4:
    # lncRNA	1	11869	14409	+	NONHSAT148173.1	.	NONHSAG000001.2	1-11869-12227,2-12613-12721,3-13221-14409	NONCODE
    for i in f:
        line=i.strip('\n').split('\t')
        chrid=line[1]
        ncRNA_ID=line[5]
        ncRNA_strand=line[4]
        # 不考虑没有strand的lncRNA
        if ncRNA_strand == '.':
            continue
        exons=line[-2].split(',')
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
