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


# lists = [AFR, AMR, EUR, EAS, SAS]
# num_lists = len(lists)

# for i in range(num_lists):
#     for j in range(i+1, num_lists):
#         count = len(set(lists[i]).intersection(set(lists[j])))
#         print(f"Number of items in common between  and : {count}")

# Number of items in common between  and : 18977366
# Number of items in common between  and : 13055537
# Number of items in common between  and : 9689233
# Number of items in common between  and : 11593198
# Number of items in common between  and : 14790484
# Number of items in common between  and : 9527984
# Number of items in common between  and : 12276146
# Number of items in common between  and : 9288814
# Number of items in common between  and : 12581143
# Number of items in common between  and : 10459877

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

# 1000G unique   pop percent notUniuq 
# 37321404 14581226 AFR.snv 39% 22740178
# 27080886 3585922 AMR.snv 13% 23494964
# 22866988 6580684 SAS.snv 29% 16286304
# 21159071 3248622 EUR.snv 15% 17910449
# 18857285 6193040 EAS.snv 33% 12664245

# 手动输出

