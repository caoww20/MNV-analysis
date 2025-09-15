import os,sys

file = sys.argv[1] # 'lncRNA_RNAsnp.txt'
url = sys.argv[2] # './'
# C1651A-G1650A   rs534252|MNV01263504|NONHSAT225272.1|+ AACTTGCCGTCA
seq_path = sys.argv[2] + "tempseq"
mnv_path = sys.argv[2] + "tempmnv"
res_path = sys.argv[2] + "res_RNAsnp.txt"

mnv_l = open(file)
for line in mnv_l:
    line= line.split('\t')
    with open(seq_path, "w") as seq:
        seq.write('>' + line[1] + '\n' + line[2] + '\n')
    with open(mnv_path, "w") as mnv:
        mnv.write(line[0] + '\t' +line[1]  + '\n')
    len_seq = len(line[2])
    if len_seq<1000:
        RNAsnp_mode=1
        if 200<=len_seq<400:
            window=100
        elif 400<=len_seq<600:
            window=200
        elif 600<=len_seq<800:
            window=300
        elif 800<=len_seq<1000:
            window=400
    else:
        RNAsnp_mode=2
        if 1000<=len_seq<1200:
            window=500
        elif 1200<=len_seq<1400:
            window=600
        elif 1400<=len_seq<1600:
            window=700
        elif len_seq>=1600:
            window=800
    command = "RNAsnp -s {tempsnp} -w {window} -m {RNAsnp_mode} -f {tempseq} | tail -n1".format(tempsnp=mnv.name,\
                                                                                                tempseq=seq.name,window=window,RNAsnp_mode=RNAsnp_mode)                                                                                      
    res = os.popen(command).read().strip().split()
    nn = '\t'.join(res[-11::])
    with open(res_path, "a+") as res:
        res.write(nn + '\n')
