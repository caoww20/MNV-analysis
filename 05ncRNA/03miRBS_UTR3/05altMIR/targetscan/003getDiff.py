# python 003getDiff.py UTR3/res_ref.txt UTR3/res_MNV.txt UTR3/diff_MNV_vs_ref.txt
import sys
input = sys.argv[1] # input='res_ref.txt'
input2 = sys.argv[2] # input2='res_MNV.txt'
output_file = sys.argv[3] # output_file='diff_MNV_vs_ref.txt'
# ref
ref_dict={}
with open(input) as f:
    next(f)
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        ref_dict[id]=[i[3],i[4],i[8]]
# alt
alt_dict={}
with open(input2) as f:
    next(f)
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        alt_dict[id]=[i[3],i[4],i[8]]

# gain
with open (output_file,'w') as w:
    for i in alt_dict:
        if i not in ref_dict:
            w.write(i+'\t'+'\t'.join(alt_dict[i])+'\tgain'+'\n')
        else:
            # 判断区间是否有交集，如果没有交集的话输出
            if int(alt_dict[i][0])>int(ref_dict[i][1]) or int(alt_dict[i][1])<int(ref_dict[i][0]):
                w.write(i+'\t'+'\t'.join(alt_dict[i])+'\tgain'+'\n')
# loss
with open (output_file,'a+') as w:
    for i in ref_dict:
        if i not in alt_dict:
            w.write(i+'\t'+'\t'.join(ref_dict[i])+'\tloss'+'\n')
        else:
            # 判断区间是否有交集，如果没有交集的话输出
            if int(ref_dict[i][0])>int(alt_dict[i][1]) or int(ref_dict[i][1])<int(alt_dict[i][0]):
                w.write(i+'\t'+'\t'.join(ref_dict[i])+'\tloss'+'\n')
