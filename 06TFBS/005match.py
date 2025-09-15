#  python 001match.py snv.txt mnv.txt
import sys
snv_file=sys.argv[1] # snv_file='/data/jinww/mnv/analyse2/09TFBS/02res/test_snv'
mnv_file=sys.argv[2] # mnv_file='/data/jinww/mnv/analyse2/09TFBS/02res/test_mnv'
res_url=sys.argv[3] # res_url='/data/jinww/mnv/analyse2/09TFBS/02res/'
snv_dict={}
with open(snv_file) as f:
    next(f)
    for line in f:
        line=line.strip().split('\t')
        motif=line[0]
        mnvid=line[4]
        rsid=line[7]
        snv_tfbs=mnvid+'-'+motif+'-'+rsid
        snv_dict[snv_tfbs]=line[-2]
mnv_dict={}
with open(mnv_file) as f:
    next(f)
    for line in f:
        line=line.strip().split('\t')
        motif=line[0]
        mnvid=line[4]
        mnv_tfbs=mnvid+'-'+motif
        rsIDs=line[7].split(',')
        rsIDs=[mnvid+'-'+motif+'-'+i for i in rsIDs]
        rsIDs_flag=[]
        for i in rsIDs:
            if i in snv_dict:
                rsIDs_flag.append(snv_dict[i])
            else:
                rsIDs_flag.append('None')
        mnv_dict[mnv_tfbs]=[line[-2]]+rsIDs_flag

## Generate mapping table: SNV-TFBS => MNV-TFBS
with open(res_url+'snv_mnv_tfbs.txt','w') as f:
    for i in mnv_dict:
        mnv_type=mnv_dict[i][0]
        snv_type=mnv_dict[i][1:]
        if mnv_type in snv_type:
            mnv_type='no change'
        for ii in snv_type:
            f.write(ii+'\t'+mnv_type+'\n')
## Generate mapping table: SNV-TFBS => MNV-TFBS
## Generate MNV-TFBS summary file
with open(res_url+'mnv_tfbs.txt','w') as f:
    for i in mnv_dict:
        mnv_type=mnv_dict[i][0]
        snv_type=mnv_dict[i][1:]
        if mnv_type not in snv_type:
            f.write(mnv_type+'\tyes\t'+'|'.join(snv_type)+'\n')
        else:
            f.write(mnv_type+'\tno\t'+'|'.join(snv_type)+'\n')
## Generate MNV-TFBS summary file


