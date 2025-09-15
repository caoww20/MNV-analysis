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
def getSNVPos(seq,pos_real,pos,alt,ncRNA_strand):
    if ncRNA_strand == '+':
        seq=list(seq)
        pos=int(pos)
        pos_index=pos_real.index(pos)
        ref_codon=seq[pos_index]
        alt_condon=alt
        res=ref_codon+str(pos_index+1)+alt_condon
    else:
        seq=toTrans(seq)
        seq=list(seq)
        pos_real_rev=pos_real[:]
        pos_real_rev.reverse()
        pos=int(pos)
        pos_index=pos_real_rev.index(pos)
        ref_codon=seq[pos_index]
        alt_condon=toTrans(alt)
        res=ref_codon+str(pos_index+1)+alt_condon
    return(res)
# 根据MNV突变位点获取改变后的序列
def getMNVPos(seq,pos_real,pos_list,alt_list,ncRNA_strand):
    if ncRNA_strand == '+':
        seq=list(seq)
        pos_list=[int(i) for i in pos_list]
        res=[]
        for i in range(len(pos_list)):
            pos_index=pos_real.index(pos_list[i])
            ref_codon=seq[pos_index]
            alt_condon=alt_list[i]
            res.append(ref_codon+str(pos_index+1)+alt_condon)
        res='-'.join(res)
    else:
        seq=toTrans(seq)
        seq=list(seq)
        pos_real_rev=pos_real[:]
        pos_real_rev.reverse()
        pos_list=[int(i) for i in pos_list]
        pos_list.reverse()
        alt_list=[toTrans(i) for i in alt_list]
        alt_list.reverse()
        res=[]
        for i in range(len(pos_list)):
            pos_index=pos_real_rev.index(pos_list[i])
            ref_codon=seq[pos_index]
            alt_condon=alt_list[i]
            res.append(ref_codon+str(pos_index+1)+alt_condon)
        res='-'.join(res)
    return(res)
####################################################################
# python 01getSeq.py  res_adjust.txt  chr_fix.fa ./

input = sys.argv[1] # input='lncRNA_adjust.txt'
region_url = sys.argv[2] # region_url='/data/jinww/mnv/MNVAnno/database/human/hg38_lncRNA.txt'
fa_url = sys.argv[3] # fa_url='chr_fix.fa'
output_file = sys.argv[4] # output_file='./'

# C1651A-G1650A   rs1,rs2-mnvid-trans_id AACTTGCCGTCA
# C1651A  rs1-mnvid-trans_id AACTTGCCGTCA
# 第一列中的位置是序列的位置，格式为GTF不是bed格式，即[]
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

# Load region
ncRNA={}
with open(region_url) as f:
    # lncRNA	1	11869	14409	+	NONHSAT148173.1	.	NONHSAG000001.2	1-11869-12227,2-12613-12721,3-13221-14409	NONCODE
    for i in f:
        line=i.strip('\n').split('\t')
        ncRNA[line[5]]=line[8]


# Obtain 5'utr and cds sequences
with open(input) as f, open(output_file+'lncRNA_RNAsnp.txt','w') as result1:
    for i in f:
        line=i.strip('\n').split('\t')
        chrid=line[0]
        mnv_pos=line[1].split(',')
        mnvid=line[2]
        ref=line[3].split(',')
        alt=line[4].split(',')
        rsids=line[5].split(',')
        # 不考虑没有rsid的SNP
        if '.' in rsids:
            continue
        info=line[6].split(' ')
        ncRNA_pos=[int(info[3]),int(info[4])]
        ncRNA_ID=info[6]
        ncRNA_strand=info[5]
        exons=ncRNA[ncRNA_ID].split(',')
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
        ref=ref.upper().replace("T", "U")

        # 获取SNV的seq
        for ii in range(len(mnv_pos)):
            alt_snv=getSNVPos(seq,pos_real,mnv_pos[ii],alt[ii],ncRNA_strand)
            newa=[alt_snv.replace("T", "U"),'|'.join([rsids[ii],mnvid,ncRNA_ID,ncRNA_strand]),ref]
            result1.write('\t'.join(newa)+'\n')
        # 获取MVN的seq
        alt_mnv=getMNVPos(seq,pos_real,mnv_pos,alt,ncRNA_strand)
        newa=[alt_mnv.replace("T", "U"),'|'.join([line[5],mnvid,ncRNA_ID,ncRNA_strand]),ref]
        result1.write('\t'.join(newa)+'\n')
####################################################################
