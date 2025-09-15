import sys

file1=sys.argv[1]
file2=sys.argv[2]
resfile=sys.argv[3]
# mnvid	refs	alts	distance	mnvtype	snvid	snv	mnv	rate
# MNV10282333	T,T	C,G	3	2	10:100003391:T:C,10:100003394:T:G	1,1	1	1.0
f=open(file1)
s1=f.readlines()
s1=s1[1:]
f.close()
f=open(file2)
s2=f.readlines()
s2=s2[1:]
f.close()

mnv={}
for i in s1:
    a=i.strip().split('\t')
    mnv[a[0]]=a

for i in s2:
    a=i.strip().split('\t')
    if a[0] not in mnv:
        mnv[a[0]]=a
    else:
        b=mnv[a[0]]
        # snv mnv rate
        ac1=a[6].split(',')
        ac2=b[6].split(',')
        snv_ac=[]
        for ii in range(len(ac1)):
            snv_ac.append(str(int(ac1[ii])+int(ac2[ii])))
        snv_ac=','.join(snv_ac)  
        mnv_ac=int(a[7])+int(b[7])
        rate=mnv_ac/int(snv_ac.split(',')[0])
        newa=a[0:6]+[snv_ac,str(mnv_ac),str(rate)]
        mnv[a[0]]=newa

f=open(resfile,'w')
f.write('mnvid\trefs\talts\tdistance\tmnvtype\tsnvid\tsnv\tmnv\trate\n')
for i in mnv:
    f.write('\t'.join(mnv[i])+'\n')
