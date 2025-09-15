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

input = sys.argv[1] # res_adjust.txt
fa_url = sys.argv[2] # chr_fix.fa
output_file = sys.argv[3] # ./


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
# rs1325725857_MNV01263504_MNV001_+
with open(input) as f, open(output_file+'preMIR_ref.seq','w') as result1,open(output_file+'preMIR_SNV.seq','w') as result2, open(output_file+'preMIR_MNV.seq','w') as result3:
    for i in f:
        line=i.strip('\n').split('\t')
        chrid=line[0]
        mnv_pos=line[1].split(',')
        mnvid=line[2]
        ref=line[3].split(',')
        alt=line[4].split(',')
        # 不考虑没有rsid的情况
        if '.' in line[5]:
            continue
        rsids=line[5].split(',')
        info=line[-1].split(' ')
        pre_miRNA_pos=[int(info[3]),int(info[4])]
        pre_miRNA_strand=info[5]
        pre_miRNA_ID=info[6]
        pre_miRNA_ID2=info[7]

        # 提取序列
        seq=chr_seq[chrid][(int(pre_miRNA_pos[0])-1):int(pre_miRNA_pos[1])]
        seq=seq.upper()
        pos_real=list(range(int(pre_miRNA_pos[0]),int(pre_miRNA_pos[1])+1))

        # 获取ref的seq
        ref=seq
        if pre_miRNA_strand == '-':
            ref = toTrans(ref)
        result1.write('>'+'_'.join([line[5],mnvid,pre_miRNA_ID,pre_miRNA_strand])+'\n')
        result1.write(ref.replace("T", "U")+'\n')

        # 获取SNV的seq
        for ii in range(len(mnv_pos)):
            alt_snv=getSNVSeq(seq,pos_real,mnv_pos[ii],alt[ii])
            if pre_miRNA_strand == '-':
                alt_snv=toTrans(alt_snv)
            result2.write('>'+'_'.join([rsids[ii],mnvid,pre_miRNA_ID,pre_miRNA_strand])+'\n')
            result2.write(alt_snv.replace("T", "U")+'\n')
        # 获取MVN的seq
        alt_mnv=getMNVSeq(seq,pos_real,mnv_pos,alt)
        if pre_miRNA_strand == '-':
            alt_mnv=toTrans(alt_mnv)
        result3.write('>'+'_'.join([line[5],mnvid,pre_miRNA_ID,pre_miRNA_strand])+'\n')
        result3.write(alt_mnv.replace("T", "U")+'\n')
####################################################################
#