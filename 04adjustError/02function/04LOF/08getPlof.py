## 载入plof基因##########################################
import pandas as pd
import numpy as np
## 载入plof基因注释库
f=open('/data/jinww/04reference/publicDB/gnomAD/plof/gnomad.v2.1.1.lof_metrics.by_transcript.txt')
plof_trans=[]
plof=[]
for i in f:
    # A1BG    ENST00000263100 true
    a=i.strip().split('\t')
    if a[2]=='true' and a[30]!='NA':
        plof_trans.append(a[1])
        plof.append(float(a[30]))
f.close()
df=pd.DataFrame(data={'trans':plof_trans,'plof':plof})
df["oe_decile"] = df.plof.rank(pct=True)
df.head()

## 载入lofeff的数据
lof={}
f=open('/data/jinww/mnv/analyse3/04adjustError/02function/04LOF/results')
for i in f:
    if i[0]!='#':
        a=i.strip().split('\t')
        snvid=a[2]
        info=a[-1].split(',')
        for ii in info:
            for iii in ii.split('|'):
                if ii.split('|')[-4] =='HC':
                    if ii.split('|')[6] not in lof:
                        lof[snvid]=[ii.split('|')[6]+'|'+ii.split('|')[-4]]
                    else:
                        lof[snvid].append(ii.split('|')[6]+'|'+ii.split('|')[-4])
f.close()

def isLOF(rsid,lof):
    flag='-'
    for k in rsid:
        if k in lof:
            flag='Rescued PTV'
            break
    return(flag)


f=open("/data/jinww/mnv/analyse3/04adjustError/02function/02canonicalTrans/merge_total.txt")
s=f.readlines()
f.close()
res=open("/data/jinww/mnv/analyse3/04adjustError/02function/04LOF/plofTrans.txt",'w')
res.write('mnvid\ttrans\ttype\n')
for i in s[1:]:
    a=i.strip().split('\t')
    mnvid=a[2]
    rsid=a[5]
    trans=a[8]
    MNVAAC=a[14].split('/')[1]
    snv1ACC=a[18].split('/')[1]
    snv2ACC=a[22].split('/')[1]
    snvACC=[snv1ACC,snv2ACC]
    if a[26]!='.':
        snv3ACC=a[26].split('/')[1]
        snvACC.append(snv3ACC)
    change_type='NA'
    # 如果MNVACC为*，snvACC没有一个*，则为Gained PTV
    if MNVAAC=='*' and '*' not in snvACC:
        change_type='Gained PTV'
    # 如果MNVACC不为*，snvACC有一个*，则为Rescued PTV
    if MNVAAC!='*' and '*' in snvACC:
        change_type='Rescued PTV'
        # change_type=isLOF(a[5].split(','),lof)
    # 如果change_type为NA
    if change_type=='NA':
        # 如果MNVACC不在snvACC中，则为missense
        if MNVAAC not in snvACC:
            change_type='missense'
        # 如果MNVACC在snvACC中，则为Unchanged
        else:
            change_type='Unchanged'
    # 如果change_type不为-则输出
    if change_type!='-':
        res.write(mnvid+'\t'+trans+'\t'+change_type+'\n')
res.close()


## 进行富集分析
def get_oe_decile(trid, oe_table=df):
    l = oe_table[oe_table.trans==trid]
    if l.shape[0]==0: return (np.nan) #for flagging that it does not exist
    else: return l.oe_decile.values[0]



df2 = pd.read_csv("/data/jinww/mnv/analyse3/04adjustError/02function/04LOF/plofTrans.txt", sep="\t")
g = df2[(df2.type=="Gained PTV")]["trans"] #gained
rr = df2[(df2.type=="Rescued PTV")]["trans"]#real rescue
un = df2[df2.type=="Unchanged"]["trans"]
m = df2[df2.type.apply(lambda x: "missense" in x)]["trans"]
alls = df2["trans"]

oe_g = g.apply(lambda x: get_oe_decile(x))
oe_rr = rr.apply(lambda x: get_oe_decile(x))
oe_un = un.apply(lambda x: get_oe_decile(x))
oe_m = m.apply(lambda x: get_oe_decile(x))
oe_all = alls.apply(lambda x: get_oe_decile(x))

#fraction of falling in top 20 percent
hp_g = sum(oe_g<0.2) / len(oe_g) 
hp_rr = sum(oe_rr<0.2) / (len(oe_rr))
hp_un = sum(oe_un<0.2) / len(oe_un)
hp_m = sum(oe_m<0.2) / len(oe_m)
# hp_all = sum(oe_all<0.2) / len(oe_all)
fracs = np.array([hp_g, hp_m,hp_un, hp_rr])
ns = np.array([len(oe_g),len(oe_m),len(oe_un),len(oe_rr)]) 
errs = np.sqrt((fracs*(1-fracs)) /ns) 

#unconstrained 50%
hp_g = sum(oe_g>0.5) / len(oe_g) 
hp_rr = sum(oe_rr>0.5) / len(oe_rr)
hp_un = sum(oe_un>0.5) / len(oe_un)
hp_m = sum(oe_m>0.5) / len(oe_m)
fracs = np.array([hp_g, hp_m,hp_un, hp_rr]) #[0.44897959, 0.33403843, 0.35262192, 0.28996283] [0.42424242, 0.37116387, 0.37346867, 0.35874439]
ns = np.array([len(oe_g),len(oe_m),len(oe_un),len(oe_rr)]) #[  490, 17019, 13082,  1345] [ 132, 3454, 2857,  223]
errs = np.sqrt((fracs*(1-fracs)) /ns) #[0.02246979, 0.0036154 , 0.00417731, 0.0123723 ] [0.04301698, 0.00822035, 0.00904989, 0.03211853]

from scipy import stats
def fisher_OR_and_pval(x1, x2, y1, y2): #1 case, 1 alt, 2 case, 2 alt
    oddsratio, pvalue = stats.fisher_exact([[x1, x2], [y1, y2]])
    return (oddsratio, pvalue)

#gained, constrained genes
frac_all = sum(oe_all<0.2) / sum(~np.isnan(oe_all)) 
(x1,x2) = (sum(oe_g<0.2)*(1+frac_all), len(oe_g))
(y1,y2) = (sum(oe_all<0.2)*(1+frac_all) -x1, len(oe_all)-x2) 
(OR, pval) = fisher_OR_and_pval(x1, x2, y1, y2)
print ("gained, constrained genes, p-val: {0}".format(pval)) #p-val: 2.1310468718606426e-09 p-val: 0.00039428218609345415

#rescued, constrained genes
frac_all = sum(oe_all<0.2) / sum(~np.isnan(oe_all)) 
(x1,x2) = (sum(oe_rr<0.2)*(1+frac_all), len(oe_rr))
(y1,y2) = (sum(oe_all<0.2)*(1+frac_all) -x1, len(oe_all)-x2)
(OR, pval) = fisher_OR_and_pval(x1, x2, y1, y2)
print ("rescued, constrained genes, p-val: {0}".format(pval)) #p-val: 0.11534053862382684 p-val: 0.5559323616089077

#gained, unconstrained genes
frac_all = sum(oe_all>0.5) / sum(~np.isnan(oe_all)) 
(x1,x2) = (sum(oe_g>0.5)*(1+frac_all), len(oe_g))
(y1,y2) = (sum(oe_all>0.5)*(1+frac_all) -x1, len(oe_all)-x2) 
(OR, pval) = fisher_OR_and_pval(x1, x2, y1, y2)
print ("gained, unconstrained genes, p-val: {0}".format(pval)) #p-val: 0.0001856396618691658 0.43325118403180907

#rescued, unconstrained genes
frac_all = sum(oe_all>0.5) / sum(~np.isnan(oe_all)) 
(x1,x2) = (sum(oe_rr>0.5)*(1+frac_all), len(oe_rr))
(y1,y2) = (sum(oe_all>0.5)*(1+frac_all) -x1, len(oe_all)-x2)
(OR, pval) = fisher_OR_and_pval(x1, x2, y1, y2)
print ("rescued, unconstrained genes, p-val: {0}".format(pval)) #p-val:0.00047430115729222384 p-val: 0.6908663491435227


# gained, constrained genes, p-val: 1.4033440383174333e-11
# rescued, constrained genes, p-val: 0.03897486996043722
# gained, unconstrained genes, p-val: 7.036968775681469e-05
# rescued, unconstrained genes, p-val: 0.6788121390564654