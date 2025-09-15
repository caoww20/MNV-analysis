# python 003mergeRes.py miranda/UTR3/diff_MNV_vs_ref.txt targetscan/UTR3/diff_MNV_vs_ref.txt merge/UTR3/diff_MNV_vs_ref.txt
import sys
input1 = sys.argv[1] # input1='miranda/UTR3/diff_MNV_vs_ref.txt'
input2 = sys.argv[2] # input2='targetscan/UTR3/diff_MNV_vs_ref.txt'
output_file = sys.argv[3] # output_file='merge/UTR3/diff_MNV_vs_ref.txt'


miranda=[]
miranda_dict={}
with open(input1) as f:
    for line in f:
        line=line.strip('\n').split('\t')
        miranda.append(line[0]+'|'+line[1])
        if line[0]+'|'+line[1] not in miranda_dict:
            miranda_dict[line[0]+'|'+line[1]]=[line]
        else:
            miranda_dict[line[0]+'|'+line[1]].append(line)
miranda=set(miranda)

targetscan=[]
targetscan_dict={}
with open(input2) as f:
    for line in f:
        line=line.strip('\n').split('\t')
        targetscan.append(line[0]+'|'+line[1])
        if line[0]+'|'+line[1] not in targetscan_dict:
            targetscan_dict[line[0]+'|'+line[1]]=[line]
        else:
            targetscan_dict[line[0]+'|'+line[1]].append(line)
targetscan=set(targetscan)

# Get intersection of the two sets
intersect=miranda.intersection(targetscan)

# Output merged results
with open(output_file,'w') as result:
   for id in intersect:
        miranda_info=miranda_dict[id]
        targetscan_info=targetscan_dict[id]
        for i in targetscan_info:
            tar_start=int(i[2])
            tar_end=int(i[3])
            tar_flag=i[5]
            for j in miranda_info:
                mir_start=int(j[6])
                mir_end=int(j[7])
                mir_flag=j[-1]
                # If intervals overlap and flags match, output the merged record
                if tar_start<=mir_end and tar_end>=mir_start:
                    if tar_flag==mir_flag:
                        newa=j[:-1] + i[4:]
                        result.write('\t'.join(newa)+'\n')
                        break

            
       