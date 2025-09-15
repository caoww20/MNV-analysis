# python 003getDiff.py res_ref_adjust.txt res_MNV_adjust.txt diff_MNV_vs_ref.txt
import sys
input = sys.argv[1] # input='res_ref_adjust.txt'
input2 = sys.argv[2] # input2='res_MNV_adjust.txt'
output_file = sys.argv[3] # output_file='diff_MNV_vs_ref.txt'
# ref
ref_dict={}
ref=[]
with open(input) as f:
    next(f)
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        ref.append(id)
        ref_dict[id]=line.strip('\n')
ref=set(ref)
# alt
alt_dict={}
alt=[]
with open(input2) as f:
    next(f)
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        alt.append(id)
        alt_dict[id]=line.strip('\n')
alt=set(alt)
# gain
gain_set=alt-ref
with open (output_file,'w') as w:
    for i in gain_set:
        w.write(alt_dict[i]+'\tgain'+'\n')
# loss
loss_set=ref-alt
with open (output_file,'a+') as w:
    for i in loss_set:
        w.write(ref_dict[i]+'\tloss'+'\n')