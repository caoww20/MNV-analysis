import math,sys
file1=sys.argv[1] #snv.txt
file2=sys.argv[2] #mnv.txt
resname=sys.argv[3] #merge.txt
#  python 04merge.py snv.txt  mnv.txt merge.txt
def my_mnv_category(snp_number,snp1_con, snp2_con,snp3_con, mnv_con, aa1, aa2, aa3,aa4):
    # 'gain stop' 'complete gain stop' 'rescue stop' "complete rescue stop" 'gain missense' 'complete gain missense'
    # 'change missense' 'complete change missense' 'lose missense' 'complete lose missense'
    res=[0,0,0,0,0,0,0,0,0,0] 
    if snp_number==2:
        type=[snp1_con, snp2_con]
        # 获得stop"!stop>=1 ->stop"
        if len([idx for idx, i in enumerate(type) if i !='stop_gain'])>=1 and mnv_con=='stop_gain':
            res[0]=1
        # 完全获得stop"!stop ->stop"
        if snp1_con!='stop_gain' and snp2_con!='stop_gain' and mnv_con=='stop_gain':
            res[1]=1

        # 拯救stop"stop>=1 -> !stop"
        if len([idx for idx, i in enumerate(type) if i =='stop_gain'])>=1 and mnv_con!='stop_gain':
            res[2]=1
        # 完全拯救stop"rescued stop:!stop ->stop"
        if snp1_con=='stop_gain' and snp2_con=='stop_gain' and mnv_con!='stop_gain':
            res[3]=1

        # 获得错译突变"至少一个syn->missense"
        if len([idx for idx, i in enumerate(type) if i =='synonymous_variant'])>=1 and mnv_con=='missense_variant':
            res[4]=1
        # 完全获得错译突变"syn+syn->missense;syn+syn+syn->missense"
        if snp1_con=='synonymous_variant' and snp2_con=='synonymous_variant' and mnv_con=='missense_variant':
            res[5]=1

        # 错误突变改变"至少一个miss A->miss B"
        if 'missense_variant' in type and mnv_con=='missense_variant':
            if snp1_con =='missense_variant':
                if aa1 != aa4:
                    res[6]=1
            if snp2_con =='missense_variant':
                if aa2 != aa4:
                    res[6]=1            
        # 完全错误突变改变"Miss A+B ->Miss C;Miss A+B+C ->Miss D"
        if snp1_con=='missense_variant' and snp2_con=='missense_variant' and mnv_con=='missense_variant':
            if aa1!=aa4 and aa2!=aa4:
                res[7]=1     

        # 错误突变丢失"至少一个Miss A->syn"
        if 'missense_variant' in type and mnv_con=='synonymous_variant':
            res[8]=1
        # 完全错误丢失"Miss A+MissB ->syn"
        if snp1_con=='missense_variant' and snp2_con=='missense_variant' and mnv_con=='synonymous_variant':
            res[9]=1
    if snp_number==3:
        type=[snp1_con, snp2_con,snp3_con]
        # 获得stop"!stop>=1 ->stop"
        if len([idx for idx, i in enumerate(type) if i !='stop_gain'])>=1 and mnv_con=='stop_gain':
            res[0]=1
        # 完全获得stop"!stop ->stop"
        if snp1_con!='stop_gain' and snp2_con!='stop_gain' and snp3_con!='stop_gain' and mnv_con=='stop_gain':
            res[1]=1

        # 拯救stop"stop>=1 -> !stop"
        if len([idx for idx, i in enumerate(type) if i =='stop_gain'])>=1 and mnv_con!='stop_gain':
            res[2]=1
        # 完全拯救stop"rescued stop:!stop ->stop"
        if snp1_con=='stop_gain' and snp2_con=='stop_gain' and snp3_con=='stop_gain' and mnv_con!='stop_gain':
            res[3]=1

        # 获得错译突变"至少一个syn->missense"
        if len([idx for idx, i in enumerate(type) if i =='synonymous_variant'])>=1 and mnv_con=='missense_variant':
            res[4]=1
        # 完全获得错译突变"syn+syn->missense;syn+syn+syn->missense"
        if snp1_con=='synonymous_variant' and snp2_con=='synonymous_variant' and snp3_con=='synonymous_variant' and mnv_con=='missense_variant':
            res[5]=1

        # 错误突变改变"至少一个miss A->miss B"
        if 'missense_variant' in type and mnv_con=='missense_variant':
            if snp1_con =='missense_variant':
                if aa1 != aa4:
                    res[6]=1
            if snp2_con =='missense_variant':
                if aa2 != aa4:
                    res[6]=1 
            if snp3_con =='missense_variant':
                if aa3 != aa4:
                    res[6]=1                                
        # 完全错误突变改变"Miss A+B ->Miss C;Miss A+B+C ->Miss D"
        if snp1_con=='missense_variant' and snp2_con=='missense_variant' and snp3_con=='missense_variant' and mnv_con=='missense_variant':
            if aa1!=aa4 and aa2!=aa4 and aa3!=aa4:
                res[7]=1     

        # 错误突变丢失"至少一个Miss A->syn"
        if 'missense_variant' in type and mnv_con=='synonymous_variant':
            res[8]=1
        # 完全错误丢失"Miss A+MissB ->syn"
        if snp1_con=='missense_variant' and snp2_con=='missense_variant' and snp3_con=='missense_variant' and mnv_con=='synonymous_variant':
            res[9]=1
    return(res)


## 读入数据
f=open(file1)
mysnv={}
# 1       115683985       G       C       ENSG00000173218 VANGL1  ENST00000355485 exon6   G988C   A330P   missense_variant        gct/Cct A/P
for i in f:
    a = i.strip('\n').split('\t')
    mysnv[':'.join(a[0:4]+[a[6]])]=a
f.close()
# mnv
f=open(file2)
s=f.readlines()
f.close()
# match
f1=open(resname,'w')
# myhead=['chr','pos','ref','alt','gene','gene_symbol','trans','exon','codon','AA','MNVType','MNVCodonC','MNVAAC','MNVLof','snv1','snv1Type','snv1CodonC','snv1AAC','snv1Lof',
# 'snv2','snv2Type','snv2CodonC','snv2AAC','snv2Lof','snv3','snv3Type','snv3CodonC','snv3AAC','snv3Lof',
# 'gain stop','complete gain stop','rescue stop',"complete rescue stop",'gain missense','complete gain missense',
# 'change missense','complete change missense','lose missense','complete lose missense']
myhead=['chr','pos','MNVID','ref','alt','rsid','gene','gene_symbol','trans','exon','codon','AA','MNVType','MNVCodonC','MNVAAC','snv1','snv1Type','snv1CodonC','snv1AAC',
'snv2','snv2Type','snv2CodonC','snv2AAC','snv3','snv3Type','snv3CodonC','snv3AAC',
'gain stop','complete gain stop','rescue stop',"complete rescue stop",'gain missense','complete gain missense',
'change missense','complete change missense','lose missense','complete lose missense']
f1.write('\t'.join(myhead)+'\n')
for i in s:
    a = i.strip('\n').split('\t')
    b = a[7:]
    newa=a[0:6]+b[0:5]+[b[5].replace('@',''),b[6],'.',b[5].replace('@','')[0]+'/'+b[5].replace('@','')[-1]]
    pos= a[1].split(',')
    ref= a[3].split(',')
    alt= a[4].split(',')
    k=0 # 用于查看mnv中的snv是否能有注释结果
    for ii in range(len(pos)):
        flag=':'.join([newa[0],pos[ii],ref[ii],alt[ii],newa[8]])
        if flag in mysnv:
            newa.extend([flag]+mysnv[flag][10:])
        else:
            k=1
            break
    if k ==0 : 
        if len(pos)==2:
            newa.extend(['.','.','.','.'])
        ref=newa[17].split('/')[0]
        c1=ref[0]
        c2=ref[1]
        c3=ref[2]
        for ii in newa[10].split(','):
            pos=int(ii[1:-1])
            aa_pos = math.ceil(pos/3)
            flag=pos-1-(aa_pos-1)*3
            if flag==0:
                c1=ii[-1]
            elif flag==1:
                c2=ii[-1]
            else:
                c3=ii[-1]
        alt=c1+c2+c3
        newa[11]=ref+'/'+alt
        if len(a[1].split(','))==2:
            newa.extend(my_mnv_category(2,newa[16], newa[20], newa[24], newa[12], newa[18], newa[22],newa[26], newa[14]))
            f1.write('\t'.join((str(e) for e in newa))+'\n')
        else:
            newa.extend(my_mnv_category(3,newa[16], newa[20], newa[24], newa[12], newa[18], newa[22],newa[26], newa[14]))
            f1.write('\t'.join((str(e) for e in newa))+'\n')
f.close()
f1.close()
