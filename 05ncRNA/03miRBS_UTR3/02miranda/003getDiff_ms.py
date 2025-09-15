# python 003getDiff2.py  diff_SNV_vs_ref.txt diff_MNV_vs_ref.txt diff_MNV_vs_SNV.txt
import sys
input1=sys.argv[1] # input1='diff_SNV_vs_ref.txt'
input2=sys.argv[2] # input2='diff_MNV_vs_ref.txt'
output_file=sys.argv[3] # output_file='diff_MNV_vs_SNV.txt'

# ref
ref_dict={}
with open(input1) as f:
    for line in f:
        i=line.strip('\n').split('\t')
        id='\t'.join(i[0:2])
        ref_dict[id]=i[-1]
# alt
alt_dict={}
with open(input2) as f:
    for line in f:
        i=line.strip('\n').split('\t')
        miRNA=i[1]
        id='\t'.join(i[0:2])
        alt_dict[id]=[i[-1]]
        info=i[0].split('|')
        rsids=info[0].split(',')
        for rsid in rsids:
            newid='|'.join([rsid]+info[1:])+'\t'+miRNA
            if newid in ref_dict:
                alt_dict[id].append(ref_dict[newid])
            else:
                alt_dict[id].append('none')
            
# output
with open(output_file,'w') as f:
    for i in alt_dict:
        mnv=alt_dict[i][0]
        snv=','.join(alt_dict[i][1:])
        if mnv in snv:
            f.write(i+'\t'+mnv+'\t'+snv+'\tno_change\n')
        else:
            f.write(i+'\t'+mnv+'\t'+snv+'\tchange\n')