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

input = sys.argv[1] # input='mature_miRNA_adjust.txt'
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



# Obtain seq
with open(input) as f, open(output_file+'miR_miranda_ref.seq','w') as result1,open(output_file+'miR_tar_ref.seq','w') as result2, open(output_file+'miR_miranda_MNV.seq','w') as result3, open(output_file+'miR_tar_MNV.seq','w') as result4:
# with open(input) as f:
    for i in f:
        # i='1\t1296122,1296127\tMNV01003213\tA,G\tG,A\trs372530591,rs619608\tall miRNA 1 1296110 1296129 - MIMAT0027354 hsa-miR-6726-3p MI0022571 exon1,exon1 miRBase\n'
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
        ncRNA_ID=info[6]
        ncRNA_symbol=info[7]
        ncRNA_strand=info[5]
        miRNA_start=int(info[3])
        miRNA_end=int(info[4])

        # 只考虑种子区域
        if ncRNA_strand == '+':
            seed_region=[int(miRNA_start)+1,int(miRNA_start)+7]
            if int(mnv_pos[0]) < seed_region[0] or int(mnv_pos[0]) > seed_region[1]:
                continue
            if int(mnv_pos[1]) < seed_region[0] or int(mnv_pos[1]) > seed_region[1]:
                continue
        else:
            seed_region=[int(miRNA_end)-7,int(miRNA_end)-1]
            if int(mnv_pos[0]) < seed_region[0] or int(mnv_pos[0]) > seed_region[1]:
                continue
            if int(mnv_pos[1]) < seed_region[0] or int(mnv_pos[1]) > seed_region[1]:
                continue

        # 输出哪些MNV
        print(mnvid+'\t'+ncRNA_symbol)

        # 提取序列
        seq=chr_seq[chrid][(miRNA_start-1):miRNA_end]
        pos_real=list(range(miRNA_start,miRNA_end+1))
        seq=seq.upper()

        # 获取ref的seq
        ref=seq
        if ncRNA_strand == '-':
            ref = toTrans(ref)
        # miranda
        result1.write('>'+ncRNA_symbol+'\n')
        result1.write(ref.replace("T", "U")+'\n')
        # tar
        result2.write(ncRNA_symbol+'\t'+ref[1:8].replace("T", "U")+'\t9606\n')

        # 获取MVN的seq
        alt_mnv=getMNVSeq(seq,pos_real,mnv_pos,alt)
        if ncRNA_strand == '-':
            alt_mnv=toTrans(alt_mnv)
        # miranda
        result3.write('>'+ncRNA_symbol+'\n')
        result3.write(alt_mnv.replace("T", "U")+'\n')
        # tar
        result4.write(ncRNA_symbol+'\t'+alt_mnv[1:8].replace("T", "U")+'\t9606\n')

####################################################################
