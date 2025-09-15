# python 001quantifier.py mature_seq.list SRR.fa SRR.quant
import sys
seq_file=sys.argv[1]
fa_file=sys.argv[2]
out_file=sys.argv[3]

seq_id={}
# hsa-miR-125a-5p TCGATCG
with open(seq_file) as f:
    for line in f:
        line=line.strip().split('\t')
        seq_id[line[1].replace('U','T')]=[line[0],0]

with open(fa_file) as f:
    for line in f:
        if line.startswith('@'):
            seq=next(f).strip()
            if len(seq)<18:
                continue
            if seq in seq_id:
                seq_id[seq][1]+=1

with open(out_file,'w') as f:
    for seq in seq_id:
        f.write('\t'.join([seq_id[seq][0],seq,str(seq_id[seq][1])])+'\n')



