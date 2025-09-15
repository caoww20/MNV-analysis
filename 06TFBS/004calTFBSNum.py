# python 004calTFBSNum.py geneMNVTFBSNum_hocomoco.txt geneMNVTFBSNum_hocomoco.num
import sys
input_file = sys.argv[1] # input_file='geneMNVTFBSNum_hocomoco.txt'
output_file = sys.argv[2] # output_file='geneMNVTFBSNum_hocomoco.num'

with open(input_file) as f, open(output_file,'w') as w:
    for line in f:
        line = line.strip().split('\t')
        mnvnum=len(line[1].split(','))
        if len(line[2])!=0:
            difMNVnum=len(line[2].split(','))
        else:
            difMNVnum=0
        if len(line[3])!=0:
            diffPairNum=len(line[3].split(','))
        else:
            diffPairNum=0
        if len(line[4])!=0:
            gainnum=len(line[4].split(','))
        else:
            gainnum=0
        if len(line[5])!=0:
            lossnum=len(line[5].split(','))
        else:
            lossnum=0
        w.write(line[0]+'\t'+str(mnvnum)+'\t'+str(difMNVnum)+'\t'+str(diffPairNum)+'\t'+str(gainnum)+'\t'+str(lossnum)+'\n')