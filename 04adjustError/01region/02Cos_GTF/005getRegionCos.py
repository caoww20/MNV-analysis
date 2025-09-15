import sys
cos=sys.argv[1]
input=sys.argv[2]
output=sys.argv[3]
        
type=input.split('.')[0]
# read cos file
cos_score={}
f=open(cos)
for i in f:
    a=i.strip().split('\t')
    id=a[0]+'_'+a[1]
    cos_score[id]=a[-1]
f.close()

# read input
res=open(output,'w')
f=open(input)
for i in f:
    a=i.strip().split('\t')
    id=a[0]+'_'+a[1]
    if id in cos_score:
        res.write(type+'\t'+i.strip()+'\t'+cos_score[id]+'\n')
f.close()
res.close()

