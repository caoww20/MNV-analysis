import sys
url=sys.argv[1] # '/data/jinww/mnv/analyse2/05population/01overlap/snv/'
res=sys.argv[2] # snv/overlap.txt

AFR=[]
with open(url+'AFR.snv', 'r') as file:
    line = file.readline()
    while line:
        AFR.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = file.readline()
AMR=[]
with open(url+'AMR.snv', 'r') as file:
    line = file.readline()
    while line:
        AMR.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = file.readline()
EUR=[]
with open(url+'EUR.snv', 'r') as file:
    line = file.readline()
    while line:
        EUR.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = file.readline()
EAS=[]
with open(url+'EAS.snv', 'r') as file:
    line = file.readline()
    while line:
        EAS.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = file.readline()
SAS=[]
with open(url+'SAS.snv', 'r') as file:
    line = file.readline()
    while line:
        SAS.append(line.strip().split('\t')[1])  # 追加每行内容到 A 中，并去除行末的换行符
        line = file.readline()


AFR=set(AFR)
AMR=set(AMR)
EUR=set(EUR)
EAS=set(EAS)
SAS=set(SAS)



non = AMR | EUR | EAS | SAS
difference = len(AFR - non)
print (difference) #14581226
non = AFR | EUR | EAS | SAS
difference = len(AMR - non)
print (difference) #3585922
non = AMR | AFR | EAS | SAS
difference = len(EUR - non)
print (difference) #3248622
non = AMR | EUR | AFR | SAS
difference = len(EAS - non)
print (difference) # 6193040
non = AMR | EUR | EAS | AFR
difference = len(SAS - non)
print (difference) # 6580684

