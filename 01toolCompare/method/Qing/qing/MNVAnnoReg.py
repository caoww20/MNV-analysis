import time,os,argparse,math

########Help############################################################
parser = argparse.ArgumentParser(description='this is a description')
parser.add_argument('-i',"--infile",required=True,dest="infile",help="The absolute path of the MNV file;")
parser.add_argument('-D',"--database",required=True,dest="database",help="The absolute path of the annotation database;")
parser.add_argument('-c',"--chr",required=True,dest="chr",help="Chromosomes list in MNV file; multiple alignments must be in a comma separated list, such as 1,2,3,X,Y;")
parser.add_argument('-o',"--output",required=True,dest="output",help="output filenmae;")
args=parser.parse_args()
####################################################################

########Functions############################################################
# Load database
def readAllAnno(chr_list,annotation_file_url):
    chr_anno_dict={}
    for i in chr_list:
        chr_anno_dict[i]=[]
    f=open(annotation_file_url)
    for i in f:
        a = i[0:50].split("\t")
        if a[1] in chr_list: 
            chr_anno_dict[a[1]].append(i)        
    f.close()
    return chr_anno_dict
# Load MNV files
def readAllMNV(chr_list,mnv_file_url):
    chr_mnv_dict={}
    for i in chr_list:
        chr_mnv_dict[i]=[]
    f=open(mnv_file_url)
    for i in f:
        if i[0]!='#':
            chr_mnv_dict[i[0:20].split("\t")[0]].append(i)
    f.close()
    return chr_mnv_dict
# Load Annotation for single chromosome
def readAnnotation(chr_anno_dict,chrid):
    my_data_list = []
    for i in chr_anno_dict[chrid]:
        a = i.strip('\n').split('\t')
        my_data_list.append(a) 
    return my_data_list
# Create annotation_dict
def createAnnoDict(annotation_list):
    # piRNA   1       10533   10559   1:10533-10559:hsa-7566157       .       .       piRBase .
    # ATAC    1       9926    11428   1:9926-11428    .       .       ATACdb  .
    annotation_dict = {}
    for i in annotation_list:
        annotation_dict[i[4]]=' '.join(i)
    return annotation_dict
# Create range_dict
def createRange(annotation_list):
    range_dict = {}
    temp = [annotation_list[0][2:]]
    trans_start = int(annotation_list[0][2])
    trans_end = int(annotation_list[0][3])
    for i in range(1,len(annotation_list)+1):
        if i != len(annotation_list):
            if int(annotation_list[i][2]) <= trans_end:
                temp.append(annotation_list[i][2:])
                if int(annotation_list[i][3]) > trans_end:
                    trans_end=int(annotation_list[i][3])
            else:
                range_dict[str(trans_start)+'-'+str(trans_end)]=temp
                temp=[annotation_list[i][2:]]
                trans_start=int(annotation_list[i][2])
                trans_end=int(annotation_list[i][3])
        else: 
            range_dict[str(trans_start)+'-'+str(trans_end)]=temp  
    return range_dict   
# Create point_list
def createPoint(dict):
    mylist=[]
    for i in dict:
        mylist.append(int(i.split('-')[0]))
        mylist.append(int(i.split('-')[1]))
    return mylist
# Split data to improve speed
def shrinkRange(point_list):
    point_split_dict={}
    if len(point_list)/2 <= 100:
        point_split_dict[str(point_list[0])+'-'+str(point_list[-1])]=point_list[:]
    else:
        num=math.ceil(len(point_list)/200)
        for i in range(num):
            if (i+1)*200 > len(point_list):
                tmp=point_list[i*200:len(point_list)]
            else:
                tmp=point_list[i*200:(i+1)*200]
            point_split_dict[str(tmp[0])+'-'+str(tmp[-1])]=tmp
    return point_split_dict
# Get region
def getKey(point,point_list):
    if point in point_list:
        point_index = point_list.index(point)
        if point_index%2 == 0:
            return str(point_list[point_index])+'-'+str(point_list[point_index+1]) 
        else:
            return str(point_list[point_index-1])+'-'+str(point_list[point_index])
    else:
        point_list_new=point_list[:] 
        point_list_new.append(point)
        point_list_new.sort()
        point_index = point_list_new.index(point)
        if point_index != 0 and point_index != len(point_list_new)-1:
            left_point = point_list_new[point_index-1]
            right_point = point_list_new[point_index+1]
            left_point_index = point_list.index(left_point)
            right_point_index = point_list.index(right_point)
            if left_point_index%2 == 0 and right_point_index%2 ==1:
                return str(left_point)+'-'+str(right_point)
#  Obtain information for MNV
def getInfo(start,end,c_start,c_end,range_dict,annotation_dict):
    info = []
    d = []
    d = list(set([c_start,c_end]))
    d = list(filter(None,d))
    for ii in d:
        for iii in range_dict[ii]:
            set1 = [start,int(iii[0])]
            set2 = [end,int(iii[1])]
            # Judge if intersection
            if max(set1) <= min(set2):
                if start >= int(iii[0]) and end <= int(iii[1]) :
                    info.append('all '+annotation_dict[iii[2]])
                else:
                    info.append('part '+annotation_dict[iii[2]])    
    return info
####################################################################

print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Run MNVAnnoReg.py')
########Load parameters############################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load parameters')
mnv_file_url = args.infile 
annotation_file_url= args.database
if not os.path.exists(annotation_file_url):
    print('no databases')
    exit()
output = args.output
if '/' in mnv_file_url:
    mnv_result_url = '/'.join(mnv_file_url.split('/')[:-1])+'/'+output
else:
    mnv_result_url=output
chr_list = args.chr 
chr_list=chr_list.split(',')
####################################################################

########Data loading, processing and output################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load data')
# Load database
chr_anno_dict=readAllAnno(chr_list,annotation_file_url)
# Load MNV files
chr_mnv_dict=readAllMNV(chr_list,mnv_file_url)
# Open files
result=open(mnv_result_url,'w')
# Cycle for chr
for chrid in chr_list:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Process the chr' + chrid)
    # Extract single chr
    annotation_list = readAnnotation(chr_anno_dict,chrid)
    if len(annotation_list)==0:
        # 注释
        for i in chr_mnv_dict[chrid]:
            result.write(i.strip('\n')+'\t.\n')        
        continue

    annotation_dict = createAnnoDict(annotation_list)
    range_dict = createRange(annotation_list)
    point_list = createPoint(range_dict)
    point_split_dict=shrinkRange(point_list)
    point_split_list=createPoint(point_split_dict)
    # Annotation
    for i in chr_mnv_dict[chrid]:
        # chr pos ID ref alt rsID 
        a = i.strip('\n').split('\t')[0:6]
        start = int(a[1].split(',')[0])
        end = int(a[1].split(',')[-1])
        # Start position
        b_start=getKey(start,point_split_list)
        c_start=None
        if b_start is not None:
            c_start=getKey(start,point_split_dict[b_start])
        # End position
        b_end=getKey(end,point_split_list)
        c_end=None
        if b_end is not None:
            c_end=getKey(end,point_split_dict[b_end])
        # MNV position
        if c_start is not None or c_end is not None:
            info = getInfo(start,end,c_start,c_end,range_dict,annotation_dict)
            a.append('|'.join(info))  
        else:
            a.append('.')           
        result.write('\t'.join(a)+'\n')
# Close files
result.close()
####################################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'End MNVAnnoReg.py')
