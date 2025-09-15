import re,argparse,time

########Help############################################################
parser = argparse.ArgumentParser(description='this is a description')
parser.add_argument('-u',"--url",required=True,dest="url",help="The absolute path of need files and end with '/';")
parser.add_argument('-g',"--gene",required=True,dest="gene",help="The gene anno result;")
parser.add_argument('-n',"--non",required=True,dest="non",help="The non anno result;; Multiple alignments must be in a comma separated list, such as non.txt,lnc.txt")
parser.add_argument('-r',"--reg",dest="reg",help="The reg anno result; Multiple alignments must be in a comma separated list, such as enhancer.txt,miRBS.txt",default='0')
args=parser.parse_args()
####################################################################

########Functions############################################################
# Load data
def readFile(file_url,file):
    f=open(file_url+file)
    database_dict = {}
    for i in f:
        a = i.strip('\n').split('\t')
        if a[6] != '.':
            database_dict['-'.join(a[0:6])]=a[6]
        else:
            database_dict['-'.join(a[0:6])]=''
    f.close()
    return database_dict
# Load mutiple data
def readMultiFile(file_url,file):
    file_list=file.split(',')
    file_dict={}
    file_dict=readFile(file_url,file_list[0])
    if  len(file_list) > 1:
        for i in range(1,len(file_list)):
            f=open(file_url+file_list[i])
            for ii in f:
                a = ii.strip('\n').split('\t')
                if a[6] != '.':
                    if file_dict['-'.join(a[0:6])]!='':
                        file_dict['-'.join(a[0:6])]=file_dict['-'.join(a[0:6])]+'|'+a[6]
                    else:
                        file_dict['-'.join(a[0:6])]=a[6]
            f.close()
    return file_dict  
# Merge three types
def merge_all(gene_dict,non_dict,reg_dict):
    # Add noncoding information
    for i in non_dict:
        if non_dict[i]!='':
            if gene_dict[i] == '':
                gene_dict[i] = non_dict[i]
            else:
                gene_dict[i] = gene_dict[i] + '$' + non_dict[i]
    # Makert 'intergenic'
    for i in gene_dict:
        if gene_dict[i] == '':
            gene_dict[i] = 'intergenic'
    # Add regulatory information
    for i in reg_dict:
        if reg_dict[i]!='':
            gene_dict[i] = gene_dict[i] + '$' + reg_dict[i]
    return gene_dict
# Merge 'gene' and 'non' data
def merge_gene_non(gene_dict,non_dict):
    # Add noncoding information
    for i in non_dict:
        if non_dict[i]!='':
            if gene_dict[i] == '':
                gene_dict[i] = non_dict[i]
            else:
                gene_dict[i] = gene_dict[i] + '$' + non_dict[i]
    # Makert 'intergenic'
    for i in gene_dict:
        if gene_dict[i] == '':
            gene_dict[i] = 'intergenic'
    return gene_dict
# split by mutiple symbol
def go_split(s, symbol):
    symbol = "[" + symbol + "]+"
    result = re.split(symbol, s)
    return [x for x in result if x]
# Write data
def write_data(filename,data,flag):
    if flag==1:
        result=open(filename,'w')
        for i in data:
            result.write(i.replace("-", "\t")+'\t'+data[i]+'\n')
        result.close()
    else:
        result=open(filename,'w')
        for i in data:
            mnv_head = i.replace("-", "\t")
            tmp = go_split(data[i],'$|')
            for ii in tmp:
                result.write(mnv_head+'\t'+ii+'\n')
        result.close()
####################################################################

print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Run MNVAnnoNonMerge.py')
########Load parameters############################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load parameters')
file_url=args.url
gene_file=args.gene
non_file=args.non
reg_file=args.reg
merge = file_url+'merge.txt'
merge_oneline = file_url+'merge_oneline.txt'
####################################################################

########Data loading, processing and output################################################
gene_dict=readFile(file_url,gene_file)   
non_dict=readMultiFile(file_url,non_file)
if reg_file != '0':
    reg_dict=readMultiFile(file_url,reg_file)
    gene_dict=merge_all(gene_dict,non_dict,reg_dict)
else:
    merge_gene_non(gene_dict,non_dict)
write_data(merge,gene_dict,1)
write_data(merge_oneline,gene_dict,0)
####################################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'End MNVAnnoNonMerge.py')











