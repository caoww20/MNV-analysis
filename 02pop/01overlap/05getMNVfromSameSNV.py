

import sys
url_root=sys.argv[1] # url='/data/jinww/mnv/analyse2/05population/01overlap/'
file=sys.argv[2] # /data/jinww/mnv/analyse2/02data/02adjustOrigin/hg38_humanMNV
url=url_root + 'snv/' # snv_url='/data/jinww/mnv/analyse2/05population/01overlap/snv/'

AFR=[]
with open(url+'AFR.snv', 'r') as f:
    line = f.readline()
    while line:
        AFR.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = f.readline()
AMR=[]
with open(url+'AMR.snv', 'r') as f:
    line = f.readline()
    while line:
        AMR.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = f.readline()
EUR=[]
with open(url+'EUR.snv', 'r') as f:
    line = f.readline()
    while line:
        EUR.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = f.readline()
EAS=[]
with open(url+'EAS.snv', 'r') as f:
    line = f.readline()
    while line:
        EAS.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = f.readline()
SAS=[]
with open(url+'SAS.snv', 'r') as f:
    line = f.readline()
    while line:
        SAS.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = f.readline()

AFR=set(AFR)
AMR=set(AMR)
EUR=set(EUR)
EAS=set(EAS)
SAS=set(SAS)
common_elements = AFR.intersection(AMR, EUR, EAS, SAS)
res=open(url+'common_snv','w')
for i in common_elements:
    res.write(i+'\n')
res.close()


# 载入common snv
common_snv = []
f=open(url+'common_snv')
for i in f:
    if i[0]!='.':
        common_snv.append(i.strip())
f.close()
common_snv=set(common_snv)

# 读入hg38
mydata2={}
mydata2['AFR']=[]
mydata2['AMR']=[]
mydata2['EAS']=[]
mydata2['EUR']=[]
mydata2['SAS']=[]
f=open(file)
for i in f:
    a=i.strip().split('\t')
    mnvid=a[2]
    chrid=a[0]
    pos=a[1].split(',')
    refs=a[3].split(',')
    alts=a[4].split(',')
    # 总数 ,这里需要确保点都在这个共有数据集中
    flag=1
    for ii in range(len(pos)):
        if ':'.join([chrid,pos[ii],refs[ii],alts[ii]]) not in common_snv:
            flag=0
            break
    if flag==1:
        if ';AFR:' in a[-1]:
            mydata2['AFR'].append(mnvid)
        if ';AMR:' in a[-1]:
            mydata2['AMR'].append(mnvid)
        if ';EAS:' in a[-1]:
            mydata2['EAS'].append(mnvid)
        if ';EUR:' in a[-1]:
            mydata2['EUR'].append(mnvid)
        if ';SAS:' in a[-1]:
            mydata2['SAS'].append(mnvid)
f.close()

# 输出共有的SNV
res=open(url_root+'venn/venn_1000G_common.txt','w')
for i in mydata2:
    for ii in mydata2[i]:
        res.write(i+'\t'+ii+'\n')
res.close()
