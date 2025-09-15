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
# extract ±50bp seq from MNV position
def get50seq(mnv_pos,seq,pos_real,strand):
    # flag表示MNV起始位点在3UTR的真实位置，因为现在只提取了50bp的长度，注意，正链的flag的提取使用start，负链的flag提取使用end
    # flag2：正链，start到起始的长度，负链，end到终止的长度
    if strand == '+':
        start_pos=int(mnv_pos[0])
        start_pos_index=pos_real.index(start_pos)
        if start_pos_index-50<=0:
            left=0
        else:
            left=start_pos_index-50
        if start_pos_index+50>=len(seq)-1:
            right=len(seq)-1
        else:
            right=start_pos_index+50
        flag=left
        flag2=start_pos_index
    elif strand == '-':
        start_pos=int(mnv_pos[-1])
        start_pos_index=pos_real.index(start_pos)
        if start_pos_index-50<=0:
            left=0
        else:
            left=start_pos_index-50
        if start_pos_index+50>=len(seq)-1:
            right=len(seq)-1
            flag=0
        else:
            right=start_pos_index+50   
            flag=len(seq[(right+1):])
        flag2=len(seq)-start_pos_index-1
    seq=seq[left:(right+1)]
    pos_real=pos_real[left:(right+1)]
    return([seq,pos_real,flag,flag2])
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

input = sys.argv[1] # input='UTR3_adjust.txt'
region_url = sys.argv[2] # region_url='en_GRCh38.108_gene.anno.txt'
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
    # Transcript ID (unique),Chromosome name ,Strand ,Gene ID ,Gene symbol ,Start (order) ,End ,CDS Start ,CDS End ,Number of exons ,Starts of exons ,Ends of exons ,(Start-2000)-(End+1000) ,Upstream ,Downstream ,5'UTR ,3'UTR ,Exons ,Splices
    # ENST00000361739 MT + ENSG00000198712 MT-CO2  7586    8269    7586 8269 1 7586 8269 5586-9269 5586-7585 8270-9269 1-7586-8269
    for i in f:
        line=i.strip('\n').split('\t')
        ncRNA[line[0]]=line

# Obtain 5'utr and cds sequences
with open(input) as f, open(output_file+'tar_ref.seq','w') as result1,open(output_file+'tar_SNV.seq','w') as result2, open(output_file+'tar_MNV.seq','w') as result3, open(output_file+'miranda_ref.seq','w') as result4,open(output_file+'miranda_SNV.seq','w') as result5, open(output_file+'miranda_MNV.seq','w') as result6:
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
        info=ncRNA[info[2]]
        ncRNA_ID=info[0]
        ncRNA_strand=info[2]
        exons=info[16].split(',')
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
        # 只考虑50bp的序列，这里要注意前50和后50能不能取到
        tmp=get50seq(mnv_pos,seq,pos_real,ncRNA_strand)
        seq=tmp[0]
        pos_real=tmp[1]
        flag=str(tmp[2])
        flag2=str(tmp[3])

        # 获取ref的seq
        ref=seq
        if ncRNA_strand == '-':
            ref = toTrans(ref)
        result1.write('|'.join([line[5],mnvid,ncRNA_ID,ncRNA_strand,flag,flag2])+'\t9606\t'+ref.replace("T", "U")+'\n')
        result4.write('>'+'|'.join([line[5],mnvid,ncRNA_ID,ncRNA_strand,flag,flag2])+'\n')
        result4.write(ref.replace("T", "U")+'\n')

        # 获取SNV的seq
        for ii in range(len(mnv_pos)):
            alt_snv=getSNVSeq(seq,pos_real,mnv_pos[ii],alt[ii])
            if ncRNA_strand == '-':
                alt_snv=toTrans(alt_snv)
            result2.write('|'.join([rsids[ii],mnvid,ncRNA_ID,ncRNA_strand,flag,flag2])+'\t9606\t'+alt_snv.replace("T", "U")+'\n')
            result5.write('>'+'|'.join([rsids[ii],mnvid,ncRNA_ID,ncRNA_strand,flag,flag2])+'\n')
            result5.write(alt_snv.replace("T", "U")+'\n')

        # 获取MVN的seq
        alt_mnv=getMNVSeq(seq,pos_real,mnv_pos,alt)
        if ncRNA_strand == '-':
            alt_mnv=toTrans(alt_mnv)
        result3.write('|'.join([line[5],mnvid,ncRNA_ID,ncRNA_strand,flag,flag2])+'\t9606\t'+alt_mnv.replace("T", "U")+'\n')
        result6.write('>'+'|'.join([line[5],mnvid,ncRNA_ID,ncRNA_strand,flag,flag2])+'\n')
        result6.write(alt_mnv.replace("T", "U")+'\n')
####################################################################
