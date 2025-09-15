import sys
filename=sys.argv[1]

# planA for 8bp
def judgeRepeatA(ref,alt):
    pairBase=['AA','AT','AC','AG','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']    
    singleBase=['AAAA','TTTT','CCCC','GGGG']
    ref_num=[]
    alt_num=[]
    for ii in pairBase:
        ref_num.append(len(ref.split(ii))-1)
        alt_num.append(len(alt.split(ii))-1)
    ref_num=max(ref_num)
    alt_num=max(alt_num)

    # 对于临近2
    ref_num2=[]
    alt_num2=[]
    for ii in singleBase:
        ref_num2.append(len(ref.split(ii))-1)
        alt_num2.append(len(alt.split(ii))-1)
    ref_num2=max(ref_num2)
    alt_num2=max(alt_num2)

    if ref_num2>0:
        ref_num=ref_num+1
    if alt_num>0:
        alt_num=alt_num+1

    if (ref_num>2 and alt_num>2) or (ref_num>3) or (alt_num>3):
        return 1
    else:
        return 0
    

    if (ref_num>2):
        return 1
    else:
        return 0

# planB for 10bp
def judgeRepeatB(ref_seq,alt_seq):
    pairBase=['AA','AT','AC','AG','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']
    singleBase=['AAAAAA','TTTTTT','CCCCCC','GGGGGG']

    ref_num=[]
    alt_num=[]

    ref_num=[]
    for ii in pairBase:
        ref_num.append(len(ref_seq.split(ii))-1)
        alt_num.append(len(alt_seq.split(ii))-1)
    ref_num=max(ref_num)
    alt_num=max(alt_num)

    ref_num2=[]
    alt_num2=[]
    for ii in singleBase:
        ref_num2.append(len(ref_seq.split(ii))-1)
        alt_num2.append(len(alt_seq.split(ii))-1)
    ref_num2=max(ref_num2)
    alt_num2=max(alt_num2)

    if ref_num2>0:
        ref_num=ref_num+1
    if alt_num>0:
        alt_num=alt_num+1

    if (ref_num>3 and alt_num>3) or (ref_num>4) or (alt_num>4):
    # if (ref_num>3 and alt_num>3):
        return 1
    else:
        return 0


def judgeOneStep(isrepeat,snv,rate):
    if isrepeat==1:
        return 0
    else:
        if float(rate)<0.9:
            return 0
        else:
            a=snv.split(',')
            a=[ int(x) for x in a ]
            # 判断a的差距是多少
            if (max(a)-min(a))/min(a)>0.1:
                return 0
            else:
                return 1
def judgeSNV(isrepeat,isonestep):
    if isrepeat ==1 or isonestep ==1:
        return 0
    else:
        return 1




# mnvid   refs    alts    distance        mnvtype snvid   snv     mnv     rate    ref_seq alt_seq
# MNV10282317     A,A,A   C,C,C   7       3       10:100000867:A:C,10:100000868:A:C,10:100000874:A:C      45,502,501      1       0.022222222222222223    AAAAAAAAAA      AACCAAAAAC
# mnvid   refs    alts    distance        mnvtype snvid   snv     mnv     rate    ref_seq alt_seq   snv_event one_step repeat

f=open(filename)
s=f.readlines()
f.close()
s=s[1:]

f1=open(filename.split('.')[0]+'_adjust.txt','w')
f1.write('mnvid\trefs\talts\tdistance\tmnvtype\tsnvid\tsnv\tmnv\trate\tref_seq\talt_seq\tsnv_event\tone_step\trepeat\n')
for i in s:
    a=i.strip().split('\t')
    ref=a[-2]
    alt=a[-1]
    # isRepeat=judgeRepeatA(ref,alt)
    isRepeat=judgeRepeatB(ref,alt)
    isOneStep=judgeOneStep(isRepeat,a[6],a[8])
    isSNV=judgeSNV(isRepeat,isOneStep)
    newa=a+[isSNV,isOneStep,isRepeat]
    newa = [ str(x) for x in newa ]
    f1.write('\t'.join(newa)+'\n')
f1.close()
del s
