import argparse,time

########Help############################################################
parser = argparse.ArgumentParser(description='this is a description')
parser.add_argument('-i',"--infile",required=True,dest="infile",help="The absolute path of the file, such as /home/username/mnv.txt;")
args=parser.parse_args()
####################################################################

print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Run toSNV.py')
########Load parameters############################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load parameters')
mnv_file_url = args.infile
output = mnv_file_url.replace('.txt','') + '.snv.txt'
####################################################################

########Process data################################################
f=open(mnv_file_url)
my_data = []
for i in f:
    a = i.strip('\n').split('\t')[0:6]
    chrid = a[0]
    point = a[1].split(',')
    ref = a[3].split(',')
    alt = a[4].split(',')
    rsid= a[5].split(',')
    for ii in range(len(point)):
        my_data.append('-'.join([chrid,point[ii],rsid[ii],ref[ii],alt[ii]]))
f.close()
# Unique and sort
my_data=list(set(my_data))
my_data.sort(key=lambda x:(x.split('-')[0],int(x.split('-')[1])))
# Write data
result=open(output,'w') 
result.write('#CHROM\tPOS\tID\tREF\tALT\n')
for i in my_data:
    result.write(i.replace('-','\t')+'\n')
result.close()
####################################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'End toSNV.py')