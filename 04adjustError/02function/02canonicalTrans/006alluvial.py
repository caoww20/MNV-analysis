
import sys
filename=sys.argv[1]
resname=sys.argv[2] # alluvial
## 绘制转变的流体##########################################
f=open(filename)
f1=open(resname,'w')
s=f.readlines()
s=s[1:]
for i in s:
    a = i.strip('\n').split('\t')
    flag=len(a[1].split(','))
    mnv_type=a[12]
    snv1_type=a[16]
    snv2_type=a[20]
    mnv_aa=a[14][-1]
    snv1_aa=a[18][-1]
    snv2_aa=a[22][-1]
    if flag==2:
        if mnv_type=='missense_variant':
            if mnv_aa!=snv1_aa and mnv_aa!=snv2_aa:
                f1.write('\t'.join([snv1_type,mnv_type])+'\n')
                f1.write('\t'.join([snv2_type,mnv_type])+'\n')
            else:
                if snv1_type == 'missense_variant':
                    f1.write('\t'.join([snv1_type,'part_missense_variant'])+'\n')
                else:
                    f1.write('\t'.join([snv1_type,mnv_type])+'\n')
                if snv2_type == 'missense_variant':
                    f1.write('\t'.join([snv2_type,'part_missense_variant'])+'\n')
                else:
                    f1.write('\t'.join([snv2_type,mnv_type])+'\n')
        else:
            f1.write('\t'.join([snv1_type,mnv_type])+'\n')
            f1.write('\t'.join([snv2_type,mnv_type])+'\n')
    else:
        snv3_type=a[24]
        snv3_aa=a[26][-1]
        if mnv_type=='missense_variant':
            if mnv_aa!=snv1_aa and mnv_aa!=snv2_aa and mnv_aa!=snv3_aa:
                f1.write('\t'.join([snv1_type,mnv_type])+'\n')
                f1.write('\t'.join([snv2_type,mnv_type])+'\n')
                f1.write('\t'.join([snv1_type,mnv_type])+'\n')     
            else:
                if snv1_type == 'missense_variant':
                    f1.write('\t'.join([snv1_type,'part_missense_variant'])+'\n')
                else:
                    f1.write('\t'.join([snv1_type,mnv_type])+'\n')
                if snv2_type == 'missense_variant':
                    f1.write('\t'.join([snv2_type,'part_missense_variant'])+'\n')
                else:
                    f1.write('\t'.join([snv2_type,mnv_type])+'\n')
                if snv3_type == 'missense_variant':
                    f1.write('\t'.join([snv3_type,'part_missense_variant'])+'\n')
                else:
                    f1.write('\t'.join([snv3_type,mnv_type])+'\n')     
        else:
            f1.write('\t'.join([snv1_type,mnv_type])+'\n')
            f1.write('\t'.join([snv2_type,mnv_type])+'\n')
            f1.write('\t'.join([snv1_type,mnv_type])+'\n')
f.close()
f1.close()

## 对于synonymous_variant->stop_lost 归到stop_lost->stop_lost
f=open(resname)
s=f.readlines()
f.close()
res=open(resname,'w')
for i in s:
    a=i.strip().split('\t')
    if  a[0]=='synonymous_variant' and a[1]=='stop_lost':
        res.write('stop_lost\tstop_lost\n')
    else:
        res.write(i)
res.close()