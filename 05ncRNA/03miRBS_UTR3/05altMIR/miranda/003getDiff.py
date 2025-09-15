# python 003getDiff.py res_ref_adjust.txt res_MNV_adjust.txt diff_MNV_vs_ref.txt
import sys
input = sys.argv[1] # input='res_ref_adjust.txt'
input2 = sys.argv[2] # input2='res_MNV_adjust.txt'
output_file = sys.argv[3] # output_file='diff_MNV_vs_ref.txt'
# ref
ref_dict={}
ref_dict_bak={}
with open(input) as f:
    next(f)
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        ref_dict[id]=[int(i[6]),int(i[7])]
        ref_dict_bak[id]=line.strip('\n')

# alt
alt_dict={}
alt_dict_bak={}
with open(input2) as f:
    next(f)
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        alt_dict[id]=[int(i[6]),int(i[7])]
        alt_dict_bak[id]=line.strip('\n')

# gain
with open (output_file,'w') as w:
    for i in alt_dict:
        if i not in ref_dict:
            w.write(alt_dict_bak[i]+'\tgain'+'\n')
        else:
            # 判断区间是否有交集，如果没有交集的话输出
            if alt_dict[i][0]>ref_dict[i][1] or alt_dict[i][1]<ref_dict[i][0]:
                w.write(alt_dict_bak[i]+'\tgain'+'\n')
# loss
with open (output_file,'a+') as w:
    for i in ref_dict:
        if i not in alt_dict:
            w.write(ref_dict_bak[i]+'\tloss'+'\n')
        else:
            # 判断区间是否有交集，如果没有交集的话输出
            if ref_dict[i][0]>alt_dict[i][1] or ref_dict[i][1]<alt_dict[i][0]:
                w.write(ref_dict_bak[i]+'\tloss'+'\n')

