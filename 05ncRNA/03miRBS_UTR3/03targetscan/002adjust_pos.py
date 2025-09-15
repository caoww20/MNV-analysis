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
    next(f)
    for line in f:
        line=line.strip('\n').split('\t')
        id=line[0].split('|')
        distance=int(id[-2])
        pos=int(id[-1])
        id='|'.join(id[:-3])
        start=int(line[3])
        end=int(line[4])
        if distance >50:
            # 判断[start,end]是否与[50，60]有交集
            if is_interval_intersect(start,end,50,60) == True :
                # 有交集
                start=str(int(line[3])+distance)
                end=str(int(line[4])+distance)
                result.write('\t'.join([id,line[1],start,end,line[8]])+'\n')
        else:
            # 判断[start,end]是否与[x，x+10]有交集
            if is_interval_intersect(start,end,pos,pos+10) == True :
                # 有交集
                start=str(int(line[3])+distance)
                end=str(int(line[4])+distance)
                result.write('\t'.join([id,line[1],start,end,line[8]])+'\n')