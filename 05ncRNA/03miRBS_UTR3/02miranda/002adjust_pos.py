# python 002adjust_pos.py res_ref.txt res_ref_adjust.txt
import sys
input = sys.argv[1] # input='res_ref.txt'
output_file = sys.argv[2] # output_file='res_ref_adjust.txt'
# 查看数据是否有交集
def is_interval_intersect(a, b, c, d):
    if b < c or d < a:
        return False
    return True

with open(input) as f, open(output_file,'w') as result:
    result.write('\t'.join(['target_seq','miRNA','score','energy','Q_start','Q_end','R_start','R_end','match_seq','Q_perc','R_perc','fig1','fig2','fig3'])+'\n')
    myres=[]
    flag=0
    myseq=['','','']
    for line in f:
        if 'Query:' in line:
            myseq[0]=line.strip('\n')
            flag=1
            continue
        if flag==1:
            myseq[1]=line.strip('\n')
            flag=0
            continue
        if 'Ref:' in line:
            myseq[2]=line.strip('\n')
            continue
        if '>' in line and '>>' not in line:
            tmp=line.strip('\n').replace('>','').replace(' ','\t').split('\t')+myseq
            myres.append(tmp)
            myseq=['','','']
            continue
    for line in myres:
        id=line[1].split('|')
        distance=int(id[-2])
        pos=int(id[-1])
        id='|'.join(id[:-3])
        start=int(line[6])
        end=int(line[7])
        # 调整顺序并替换id
        line=[id,line[0]]+line[2:]
        if distance >50:
            # 判断[start,end]是否与[50，60]有交集
            if is_interval_intersect(start,end,50,60) == True :
                # 有交集
                line[6]=str(start+distance)
                line[7]=str(end+distance)
                result.write('\t'.join(line)+'\n')
        else:
            # 判断[start,end]是否与[x，x+10]有交集
            if is_interval_intersect(start,end,pos,pos+10) == True :
                # 有交集
                line[6]=str(start+distance)
                line[7]=str(end+distance)
                result.write('\t'.join(line)+'\n')