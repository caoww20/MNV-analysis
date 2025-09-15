# 两个方法均为一致的记作一致,如果是一个none一个gain则为gain，如果既有loss又有gain则标记为unclear
# python 001mergeRes.py ../02miranda/UTR3/diff_MNV_vs_ref.txt ../03targetscan/UTR3/diff_MNV_vs_ref.txt UTR3/diff_MNV_vs_ref.txt &
import sys
input1 = sys.argv[1] # input1='../02miranda/UTR3/diff_MNV_vs_ref.txt'
input2 = sys.argv[2] # input2='../03targetscan/UTR3/diff_MNV_vs_ref.txt'
output_file = sys.argv[3] # output_file='UTR3/diff_MNV_vs_ref.txt'

miranda=[]
miranda_dict={}
with open(input1) as f:
    for line in f:
        line=line.strip('\n').split('\t')
        miranda.append(line[0]+'|'+line[1])
        miranda_dict[line[0]+'|'+line[1]]=line
miranda=set(miranda)

targetscan=[]
targetscan_dict={}
with open(input2) as f:
    for line in f:
        line=line.strip('\n').split('\t')
        targetscan.append(line[0]+'|'+line[1])
        targetscan_dict[line[0]+'|'+line[1]]=line[4:]
targetscan=set(targetscan)

# 获取两者的交集
intersect=miranda.intersection(targetscan)

# 获取差集
only_m = miranda - targetscan
only_t = targetscan - miranda

# 输出
with open(output_file,'w') as result:
    for id in only_m:
        miranda_info=miranda_dict[id]
        new_line=miranda_info[:-1]+['-',miranda_info[-1]]
        result.write('\t'.join(new_line)+'\n')
    for id in only_t:
        targetscan_info=targetscan_dict[id]
        new_line= ['|'.join(id.split('|')[:-1]),id.split('|')[-1]] + ['-'] * 12+ targetscan_info
        result.write('\t'.join(new_line)+'\n')
    for id in intersect:
        miranda_info=miranda_dict[id]
        miranda_flag=miranda_info[-1]
        targetscan_info=targetscan_dict[id]
        targetscan_flag=targetscan_info[-1]
        if miranda_flag == targetscan_flag:
            flag=miranda_flag
        else:
            flag='unclear'
        new_line=miranda_info[:-1]+[targetscan_info[0],flag]
        result.write('\t'.join(new_line)+'\n')