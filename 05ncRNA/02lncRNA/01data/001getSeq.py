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
# python 01getSeq.py  res_adjust.txt  chr_fix.fa ./

input = sys.argv[1] # input='res_adjust.txt'
region_url = sys.argv[2] # region_url='MNVAnno/database/human/lncRNA.txt'
fa_url = sys.argv[3] # fa_url='chr_fix.fa'
output_file = sys.argv[4] # output_file='./'


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
with open(input) as f, open(output_file+'lncRNA_ref.seq','w') as result1,open(output_file+'lncRNA_SNV.seq','w') as result2, open(output_file+'lncRNA_MNV.seq','w') as result3:
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
        result1.write('>'+'|'.join([line[5],mnvid,ncRNA_ID])+'\n')
        result1.write(ref.replace("T", "U")+'\n')

        # 获取SNV的seq
        for ii in range(len(mnv_pos)):
            alt_snv=getSNVSeq(seq,pos_real,mnv_pos[ii],alt[ii])
            if ncRNA_strand == '-':
                alt_snv=toTrans(alt_snv)
            result2.write('>'+'|'.join([rsids[ii],mnvid,ncRNA_ID])+'\n')
            result2.write(alt_snv.replace("T", "U")+'\n')

        # 获取MVN的seq
        alt_mnv=getMNVSeq(seq,pos_real,mnv_pos,alt)
        if ncRNA_strand == '-':
            alt_mnv=toTrans(alt_mnv)
        result3.write('>'+'|'.join([line[5],mnvid,ncRNA_ID])+'\n')
        result3.write(alt_mnv.replace("T", "U")+'\n')
####################################################################
