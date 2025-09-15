import sys
cos=sys.argv[1]
chrid=sys.argv[2]
input=sys.argv[3]
output=sys.argv[4]
        
# read cos file
cos_score={}
f=open(cos)
for i in f:
    a=i.strip().split('\t')
    cos_score[a[0]]=a[1]
f.close()

# read input
res=open(output,'w')
f=open(input)
for i in f:
    a=i.strip().split('\t')
    if a[0]==chrid:
        if a[1] in cos_score:
            res.write(i.strip()+'\t'+cos_score[a[1]]+'\n')
        else:
            res.write(i.strip()+'\t'+'NA'+'\n')
f.close()
res.close()

