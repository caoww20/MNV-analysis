import math,time,os,argparse

########Help############################################################
parser = argparse.ArgumentParser(description='this is a description')
parser.add_argument('-i',"--infile",required=True,dest="infile",help="The absolute path of the MNV file;")
parser.add_argument('-D',"--databse",required=True,dest="database",help="The absolute path of the annotation database;")
parser.add_argument('-c',"--chr",required=True,dest="chr",help="Chromosomes list in MNV file; Multiple alignments must be in a comma separated list, such as 1,2,3,X,Y;")
parser.add_argument('-u',"--upstream",dest="upstream",help="Range of upstream is 500~5000; Default is 2000;",default=2000)
parser.add_argument('-d',"--downstream",dest="downstream",help="Range of downstream is 500~5000; Default is 1000;",default=1000)
parser.add_argument('-s',"--splice",dest="splice",help="Range of splicing junction is 2~5; Default is 2;",default=2)
parser.add_argument('-o',"--output",required=True,dest="output",help="output filenmae;")
args=parser.parse_args()
####################################################################

########Functions############################################################
# Synchronous sorting
def sortByIndex(arr,arr_name):
    point_index=sorted(range(len(arr)), key=lambda k: arr[k])
    tmp=[]
    for i in point_index:
        tmp.append(arr_name[i])
    arr2=sorted(arr)
    return [arr2,tmp]
# Codons and amino acids
def createCodon():
    codon2aa_dict={}
    codon2aa_dict["TTT"]="F"
    codon2aa_dict["TTC"]="F"
    codon2aa_dict["TCT"]="S"
    codon2aa_dict["TCC"]="S"
    codon2aa_dict["TAT"]="Y"
    codon2aa_dict["TAC"]="Y"
    codon2aa_dict["TGT"]="C"
    codon2aa_dict["TGC"]="C"
    codon2aa_dict["TTA"]="L"
    codon2aa_dict["TCA"]="S"
    codon2aa_dict["TAA"]="*"
    codon2aa_dict["TGA"]="*"
    codon2aa_dict["TTG"]="L"
    codon2aa_dict["TCG"]="S"
    codon2aa_dict["TAG"]="*"
    codon2aa_dict["TGG"]="W"
    codon2aa_dict["CTT"]="L"
    codon2aa_dict["CTC"]="L"
    codon2aa_dict["CCT"]="P"
    codon2aa_dict["CCC"]="P"
    codon2aa_dict["CAT"]="H"
    codon2aa_dict["CAC"]="H"
    codon2aa_dict["CGT"]="R"
    codon2aa_dict["CGC"]="R"
    codon2aa_dict["CTA"]="L"
    codon2aa_dict["CTG"]="L"
    codon2aa_dict["CCA"]="P"
    codon2aa_dict["CCG"]="P"
    codon2aa_dict["CAA"]="Q"
    codon2aa_dict["CAG"]="Q"
    codon2aa_dict["CGA"]="R"
    codon2aa_dict["CGG"]="R"
    codon2aa_dict["ATT"]="I"
    codon2aa_dict["ATC"]="I"
    codon2aa_dict["ACT"]="T"
    codon2aa_dict["ACC"]="T"
    codon2aa_dict["AAT"]="N"
    codon2aa_dict["AAC"]="N"
    codon2aa_dict["AGT"]="S"
    codon2aa_dict["AGC"]="S"
    codon2aa_dict["ATA"]="I"
    codon2aa_dict["ACA"]="T"
    codon2aa_dict["AAA"]="K"
    codon2aa_dict["AGA"]="R"
    codon2aa_dict["ATG"]="M"
    codon2aa_dict["ACG"]="T"
    codon2aa_dict["AAG"]="K"
    codon2aa_dict["AGG"]="R"
    codon2aa_dict["GTT"]="V"
    codon2aa_dict["GTC"]="V"
    codon2aa_dict["GCT"]="A"
    codon2aa_dict["GCC"]="A"
    codon2aa_dict["GAT"]="D"
    codon2aa_dict["GAC"]="D"
    codon2aa_dict["GGT"]="G"
    codon2aa_dict["GGC"]="G"
    codon2aa_dict["GTA"]="V"
    codon2aa_dict["GTG"]="V"
    codon2aa_dict["GCA"]="A"
    codon2aa_dict["GCG"]="A"
    codon2aa_dict["GAA"]="E"
    codon2aa_dict["GAG"]="E"
    codon2aa_dict["GGA"]="G"
    codon2aa_dict["GGG"]="G"
    return codon2aa_dict
# Types
def createTypes():
    types_dict={}
    types_dict['up']='upstream_variant'
    types_dict['down']='downstream_variant'
    types_dict['spliceAcceptor']='splice_acceptor_variant'
    types_dict['spliceRegion']='splice_region_variant'
    types_dict['spliceDonor']='splice_donor_variant'
    types_dict['intron']='intron_variant'
    # types_dict['UTR5']='UTR5_variant'
    types_dict['UTR3']='UTR3_variant'
    types_dict['other']='other_variant'
    return types_dict
# Codon change
def createCodonChange():
    codon_change={}
    codon_change['A']='T'
    codon_change['T']='A'
    codon_change['C']='G'
    codon_change['G']='C'
    return codon_change
########read-in and adjust data#########################
# Load databases
def readAllAnno(chr_list,annotation_file_url):
    chr_anno_dict={}
    for i in chr_list:
        chr_anno_dict[i]=[]
    f=open(annotation_file_url)
    for i in f:
        if i[0:50].split("\t")[1] in chr_list:
            chr_anno_dict[i[0:50].split("\t")[1]].append(i)
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
    my_data_list = {}
    for i in chr_anno_dict[chrid]:
        a = i.strip('\n').split('\t')
        my_data_list[a[0]]=a
    return my_data_list
# adjust annotation_dict according 'upstream', 'downstream' and 'splicing junction' parameters
def fixAnnotation(upstream,downstream,splice,annotation_dict):
    up_diff=upstream-2000
    down_diff=downstream-1000
    splice_diff=splice-2
    if up_diff != 0 or down_diff !=0:
        for i in annotation_dict:
            if annotation_dict[i][2]=='+':
                tmp = annotation_dict[i][12].split('-')
                annotation_dict[i][12]= str(int(tmp[0])-up_diff)+'-'+str(int(tmp[1])+down_diff)  #total
                tmp = annotation_dict[i][13].split('-')
                annotation_dict[i][13]= str(int(tmp[0])-up_diff)+'-'+tmp[1] #up
                tmp = annotation_dict[i][14].split('-')
                annotation_dict[i][14]= tmp[0]+'-'+str(int(tmp[1])+down_diff) #down
            else:
                tmp = annotation_dict[i][12].split('-')
                annotation_dict[i][12]= str(int(tmp[0])+up_diff)+'-'+str(int(tmp[1])-down_diff)  #total
                tmp = annotation_dict[i][13].split('-')
                annotation_dict[i][13]= str(int(tmp[0])+up_diff)+'-'+tmp[1] #up
                tmp = annotation_dict[i][14].split('-')
                annotation_dict[i][14]= tmp[0]+'-'+str(int(tmp[1])-down_diff) #down      
    if splice_diff != 0:
        for i in annotation_dict:
            if annotation_dict[i][18]: #splice
                tmp = annotation_dict[i][18].split(',')
                tmp2=[]
                if annotation_dict[i][2]=='+':
                    for ii in range(len(tmp)):
                        tmp3=tmp[ii].split('-')
                        if ii == 0:
                            tmp2.append(tmp3[0]+'-'+tmp3[1]+'-'+str(int(tmp3[2])+splice_diff))
                        elif ii == len(tmp) -1:
                            tmp2.append(tmp3[0]+'-'+str(int(tmp3[1])-splice_diff)+'-'+tmp3[2])
                        else:
                            tmp2.append(str(int(tmp3[0])-splice_diff)+'-'+tmp3[1]+'-'+tmp3[2]+'-'+tmp3[3]+'-'+str(int(tmp3[4])+splice_diff))
                    annotation_dict[i][18]=','.join(tmp2) #down
                else:
                    for ii in range(len(tmp)):
                        tmp3=tmp[ii].split('-')
                        if ii == 0:
                            tmp2.append(tmp3[0]+'-'+tmp3[1]+'-'+str(int(tmp3[2])-splice_diff))
                        elif ii == len(tmp) -1:
                            tmp2.append(tmp3[0]+'-'+str(int(tmp3[1])+splice_diff)+'-'+tmp3[2])
                        else:
                            tmp2.append(str(int(tmp3[0])+splice_diff)+'-'+tmp3[1]+'-'+tmp3[2]+'-'+tmp3[3]+'-'+str(int(tmp3[4])-splice_diff))
                    annotation_dict[i][18]=','.join(tmp2) #down
    return annotation_dict
# Creat point_dict for each transcript
def createPoint(annotation_dict):
    point_dict={}
    for i in annotation_dict:
        point=[]
        name=[]
        # up
        tmp=annotation_dict[i][13].split('-')
        point.append(int(tmp[0]))
        point.append(int(tmp[1]))
        name.append('up')
        name.append('up')
        #down
        tmp=annotation_dict[i][14].split('-')
        point.append(int(tmp[0]))
        point.append(int(tmp[1]))
        name.append('down')
        name.append('down')
        #UTR5
        if annotation_dict[i][15]:
            tmp=annotation_dict[i][15].split(',')
            for ii in tmp:
                tmp2=ii.split('-')
                point.append(int(tmp2[1]))
                point.append(int(tmp2[2]))
                name.append('UTR5')
                name.append('UTR5')  
        #UTR3
        if annotation_dict[i][16]:
            tmp=annotation_dict[i][16].split(',')
            for ii in tmp:
                tmp2=ii.split('-')
                point.append(int(tmp2[1]))
                point.append(int(tmp2[2]))
                name.append('UTR3')
                name.append('UTR3')
        #exon
        tmp=annotation_dict[i][17].split(',')
        for ii in tmp:
            tmp2=ii.split('-')
            point.append(int(tmp2[1]))
            point.append(int(tmp2[2]))
            name.append('exon'+tmp2[0])
            name.append('exon'+tmp2[0])    
        #splice
        if annotation_dict[i][18]:
            tmp=annotation_dict[i][18].split(',')
            for ii in tmp:
                tmp2=ii.split('-')
                if len(tmp2)==3:
                    point.append(int(tmp2[1]))
                    point.append(int(tmp2[2]))
                    name.append('splice')
                    name.append('splice')
                else:
                    point.append(int(tmp2[0]))
                    point.append(int(tmp2[1]))
                    point.append(int(tmp2[3]))
                    point.append(int(tmp2[4]))
                    name.append('splice')
                    name.append('splice')
                    name.append('splice')
                    name.append('splice')
        point_dict[i]=sortByIndex(point,name)
    return point_dict
# Create range_dict
def createRange(annotation_dict):
    point=[]
    name=[]
    for i in annotation_dict:
        name.append(annotation_dict[i][0])
        if annotation_dict[i][2] == '+':
            point.append(int(annotation_dict[i][12].split('-')[0]))
        else:
            point.append(int(annotation_dict[i][12].split('-')[1]))
    range_list = sortByIndex(point,name)
    range_list = range_list[1]

    range_dict = {}
    trans_name = range_list[0]
    temp = [trans_name]
    if annotation_dict[trans_name][2] == '+':
        trans_start = int(annotation_dict[trans_name][12].split('-')[0])
        trans_end = int(annotation_dict[trans_name][12].split('-')[1])
    else:
        trans_start = int(annotation_dict[trans_name][12].split('-')[1])
        trans_end = int(annotation_dict[trans_name][12].split('-')[0])
    for i in range(1,len(range_list)+1):
        if i != len(range_list):
            trans_name = range_list[i]
            if annotation_dict[trans_name][2] == '+':
                a = int(annotation_dict[trans_name][12].split('-')[0])
                b = int(annotation_dict[trans_name][12].split('-')[1])
            else:
                a = int(annotation_dict[trans_name][12].split('-')[1])
                b = int(annotation_dict[trans_name][12].split('-')[0])         
            if a <= trans_end:
                temp.append(trans_name)
                if b > trans_end:
                    trans_end=b
            else:
                range_dict[str(trans_start)+'-'+str(trans_end)]=temp
                temp=[trans_name]
                trans_start=a
                trans_end=b
        else: 
            range_dict[str(trans_start)+'-'+str(trans_end)]=temp  
    return range_dict
# Create cds_dict
def createCDS(annotation_dict):
    cds_dict = {}
    for i in annotation_dict:
        cds_dict[i] = ','.join(annotation_dict[i][7:9])+','+','.join(annotation_dict[i][10:12])
    return cds_dict
# load MNV of single chromosome
def readMNV(chr_mnv_dict,chrid):
    my_data_list = []
    for i in chr_mnv_dict[chrid]:
        a = i.strip('\n').split('\t')
        my_data_list.append(a[0:6])
    return my_data_list
# Create mnv_dict
def createMNV(mnv_list):
    mnv_dict = {}
    for i in mnv_list:
        pos = i[1].split(',')
        for ii in pos:
            mnv_dict[int(ii)]=''
    return mnv_dict
########Annotation for each point#########################
# Obtain transcripts having the point
def getKey(one_point,point_list):
    one_point=int(one_point)
    if one_point in point_list:
        point_index = point_list.index(one_point)
        if point_index%2 == 0:
            return str(point_list[point_index])+'-'+str(point_list[point_index+1]) 
        else:
            return str(point_list[point_index-1])+'-'+str(point_list[point_index])
    else:
        point_list_new=point_list[:] 
        point_list_new.append(one_point)
        point_list_new.sort()
        point_index = point_list_new.index(one_point)
        if point_index != 0 and point_index != len(point_list_new)-1:
            left_point = point_list_new[point_index-1]
            right_point = point_list_new[point_index+1]
            left_point_index = point_list.index(left_point)
            right_point_index = point_list.index(right_point)
            if left_point_index%2 == 0 and right_point_index%2 ==1:
                return str(left_point)+'-'+str(right_point)
# Obtain type for the point
def getPointInfo(point,trans_name,info,strand,s,e):
    point_loc=info[0]
    point_name=info[1]
    point=int(point)
    if strand == '+':
        trans_s = int(s)
        trans_e = int(e)
    else:
        trans_s = int(e)
        trans_e = int(s)        
    if point in point_loc:
        point_index = point_loc.index(point)
        res = point_name[point_index]
        if res == 'up':
            return trans_name + '$up:'+str(abs(point-trans_s))
        if res == 'down':
            return trans_name + '$down:'+str(abs(point-trans_e))
        if res == 'splice':
            if strand =='+':
                if point_name[point_index-1] != 'splice':
                    return trans_name + '$spliceDonor'
                else:
                    if point_name[point_index-2] !='splice':
                        flag_splice=abs(point-point_loc[point_index-1])
                        if flag_splice == 1:
                            return trans_name + '$spliceDonor'
                        else:
                            return trans_name + '$spliceRegion'
                    else:
                        if point_name[point_index-3] !='splice':
                            flag_splice=abs(point-point_loc[point_index+1])
                            if flag_splice == 1:
                                return trans_name + '$spliceAcceptor'
                            else:
                                return trans_name + '$spliceRegion'
                        else:
                            return trans_name + '$spliceAcceptor'
            else:
                if point_name[point_index-1] != 'splice':
                    return trans_name + '$spliceAcceptor'
                else:
                    if point_name[point_index-2] !='splice':
                        flag_splice=abs(point-point_loc[point_index-1])
                        if flag_splice == 1:
                            return trans_name + '$spliceAcceptor'
                        else:
                            return trans_name + '$spliceRegion'
                    else:
                        if point_name[point_index-3] !='splice':
                            flag_splice=abs(point-point_loc[point_index+1])
                            if flag_splice == 1:
                                return trans_name + '$spliceDonor'
                            else:
                                return trans_name + '$spliceRegion'
                        else:
                            return trans_name + '$spliceDonor'                
        return trans_name + '$' + res
    else:
        point_loc_new=point_loc[:] 
        point_loc_new.append(point)
        point_loc_new.sort()
        point_index = point_loc_new.index(point)
        if point_index != 0 and point_index != len(point_loc_new)-1:
            left_point = point_loc_new[point_index-1]
            right_point = point_loc_new[point_index+1]
            left_point_index = point_loc.index(left_point)
            right_point_index = point_loc.index(right_point)
            if point_name[left_point_index] != point_name[right_point_index]:
                return trans_name + '$intron'
            else:
                if left_point_index%2 == 0 and right_point_index%2 ==1:
                    res = point_name[left_point_index]
                    if res == 'up':
                        return trans_name + '$up:'+str(abs(point-trans_s))
                    if res == 'down':
                        return trans_name + '$down:'+str(abs(point-trans_e))
                    if res == 'splice':
                        if strand =='+':
                            if point_name[left_point_index-1] != 'splice':
                                flag_splice=abs(point_loc[left_point_index]-point)
                                if flag_splice == 1:
                                    return trans_name + '$spliceDonor'
                                else:
                                    return trans_name + '$spliceRegion'
                            else:
                                if point_name[left_point_index-3] != 'splice':
                                    flag_splice=abs(point_loc[left_point_index+1]-point)
                                    if flag_splice ==1:
                                        return trans_name + '$spliceAcceptor'
                                    else:
                                        return trans_name + '$spliceRegion'
                        else:
                            if point_name[left_point_index-1] != 'splice':
                                flag_splice=abs(point_loc[left_point_index]-point)
                                if flag_splice == 1:
                                    return trans_name + '$spliceAcceptor'
                                else:
                                    return trans_name + '$spliceRegion'
                            else:
                                if point_name[left_point_index-3] != 'splice':
                                    flag_splice=abs(point_loc[left_point_index+1]-point)
                                    if flag_splice ==1:
                                        return trans_name + '$spliceDonor'
                                    else:
                                        return trans_name + '$spliceRegion'
                    return trans_name + '$' + res
                else:
                    return trans_name + '$intron'
# Annotate single point
def snvAnno(range_dict,mnv_dict,point_dict,annotation_dict):
    point=[]
    for i in range_dict:
        point.append(int(i.split('-')[0]))
        point.append(int(i.split('-')[1]))

    for i in mnv_dict:
        mnv_info=[] 
        a=getKey(i,point)
        if a:
            trans=range_dict[a]
            for ii in trans:
                b=getPointInfo(i,ii,point_dict[ii],annotation_dict[ii][2],annotation_dict[ii][5],annotation_dict[ii][6])
                if b:
                    mnv_info.append(b)
        if len(mnv_info) > 0:
            mnv_dict[i]=mnv_info
    return mnv_dict
########Annotation MNV#########################
# Obtain real genomic location of the codon
def getPosList(pos):
    temp=[]
    for i in range(int(len(pos)/2)):
        temp=temp+([i for i in range(pos[2*i], pos[(2*i+1)]+1, 1)])
    return temp
# Reverse list
def toReverse(lst):
    lst.reverse()
    return lst    
# AA change type
def getAAChangeType(ref_aa,alt_aa,is_start):
    if ref_aa==alt_aa:
        return 'synonymous_variant'
    else:
        if is_start == '1':
            if ref_aa == 'M' and alt_aa != 'M':
                return 'start_lost'
        if ref_aa == '*' and alt_aa != '*':
            return 'stop_lost'
        if ref_aa != '*' and alt_aa == '*':
            return 'stop_gain'
        return 'missense_variant'
# Unique and keep previous order
def getUniqAndSort(mylist):
    tmp = list(set(mylist))
    tmp.sort(key=mylist.index)
    return tmp
# Obtain base information of MNV[subname subdict trans]
def getTrans(i,mnv_dict):
    subname=[]
    pos = i[1]
    for ii in pos.split(','):
        subname.append(int(ii))
    subdict={x:mnv_dict[x] for x in subname if x in mnv_dict} 
    # extract transcripts
    trans=[]
    for value in subdict.values():
        for ii in value:
            trans.append(ii.split('$')[0])
    trans = getUniqAndSort(trans)
    return [subname,subdict,trans]
# Judge if MNV falls on exon
def judgeExonUTR5(subname,trans_name,subdict):
    flag=[0,0]
    for i in subname:
        if i in subdict:
            for ii in subdict[i]:
                if ii.split('$')[0] == trans_name :
                    if ii.split('$')[1] == 'UTR5':
                        flag[0] = 1
                    elif ii.split('$')[1][0:2] == 'ex':
                        flag[1] = 1
    return flag
# Obtain position information of transcript for MNV
def getInfo(subname,trans_name,subdict):
    temp = []
    for i in subname:
        if i in subdict:
            temp_dict = {}
            for ii in subdict[i]:
                temp_dict[ii.split('$')[0]] = ii.split('$')[1]
            if trans_name in temp_dict:
                temp.append(temp_dict[trans_name])
            else:
                temp.append('other') 
            del temp_dict  
        else:
            temp.append('other')
    return(temp)
# Judge start_gain
def getSeqInfoUTR5(i,ii,annotation_dict,position_info,subname):
    # ii means transcript name
    mnv_info = i[:]
    seq = annotation_dict[ii][19]
    utr5_region = annotation_dict[ii][15].split(',')
    pos=[]
    for iii in utr5_region:
        pos.append(iii.split('-')[1])
        pos.append(iii.split('-')[2])
    pos = list(map(int, pos))
    pos.sort()
    pos_list = getPosList(pos)
    if annotation_dict[ii][2] == '+':
        # Find region
        res_utr5 = [m for m, i in enumerate(position_info) if i == 'UTR5']
        utr5_left = subname[res_utr5[0]]
        utr5_right = subname[res_utr5[-1]]
        if utr5_left-pos[0]>=2:
            utr5_left = pos_list[pos_list.index(utr5_left)-2]
        else:
            utr5_left = pos[0]
        if pos[-1]-utr5_right>=2:
            utr5_right = pos_list[pos_list.index(utr5_right)+2]
        else:
            utr5_right=pos[-1]
        # Obtain ref seq and alt seq
        ref_seq=seq[pos_list.index(int(utr5_left)):(pos_list.index(int(utr5_right))+1)]
        alt_seq=ref_seq[:]
        pos_list=pos_list[pos_list.index(int(utr5_left)):(pos_list.index(int(utr5_right))+1)]
        for iii in range(len(position_info)):
            if position_info[iii] == 'UTR5':
                single_pos = int(subname[iii])
                single_pos_index = pos_list.index(single_pos)
                ref_seq=ref_seq[:single_pos_index] + mnv_info[3].split(',')[iii] + ref_seq[single_pos_index+1:]
                alt_seq=alt_seq[:single_pos_index] + mnv_info[4].split(',')[iii] + alt_seq[single_pos_index+1:] 
    else:
        # Convert to plus strand to deal
        mnv_info[3] = ','.join(mnv_info[3].split(',')[::-1])
        mnv_info[4] = ','.join(mnv_info[4].split(',')[::-1])
        pos = [ -x for x in pos]
        pos.sort()
        pos_list = [ -x for x in pos_list]
        pos_list.sort()
        subname_minus = [ -x for x in subname]
        subname_minus.sort()
        position_info_minus = position_info[:]
        position_info_minus.reverse()
        # Find region
        res_utr5 = [m for m, i in enumerate(position_info_minus) if i == 'UTR5']
        utr5_left = subname_minus[res_utr5[0]]
        utr5_right = subname_minus[res_utr5[-1]]
        if utr5_left-pos[0]>=2:
            utr5_left = pos_list[pos_list.index(utr5_left)-2]
        else:
            utr5_left = pos[0]
        if pos[-1]-utr5_right>=2:
            utr5_right = pos_list[pos_list.index(utr5_right)+2]
        else:
            utr5_right=pos[-1]
        # Obtain ref seq and alt seq
        ref_seq=seq[pos_list.index(int(utr5_left)):(pos_list.index(int(utr5_right))+1)]
        alt_seq=ref_seq[:]
        pos_list=pos_list[pos_list.index(int(utr5_left)):(pos_list.index(int(utr5_right))+1)]
        for iii in range(len(position_info_minus)):
            if position_info_minus[iii] == 'UTR5':
                single_pos = int(subname_minus[iii])
                single_pos_index = pos_list.index(single_pos)
                ref_seq=ref_seq[:single_pos_index] + codon_change[mnv_info[3].split(',')[iii]] + ref_seq[single_pos_index+1:]
                alt_seq=alt_seq[:single_pos_index] + codon_change[mnv_info[4].split(',')[iii]] + alt_seq[single_pos_index+1:] 
    # Judge
    ref_seq = ref_seq.replace('ATG','*')
    alt_seq = alt_seq.replace('ATG','*')
    utr5_ref_flag = [m for m, i in enumerate(ref_seq) if i == '*']
    utr5_alt_flag = [m for m, i in enumerate(alt_seq) if i == '*']
    utr5_flag='UTR5_variant'
    if len(utr5_alt_flag)-len(utr5_ref_flag) >0:
        utr5_flag = 'start_gain'
    elif len(utr5_alt_flag)-len(utr5_ref_flag) == 0:
        if utr5_alt_flag!=utr5_ref_flag:
            utr5_flag = 'start_gain'
    return utr5_flag
# Obtain results according MNV and transcript
def getSeqInfo(i,ii, annotation_dict,cds_dict,position_info,subname):
    mnv_info = i[:]
    seq = annotation_dict[ii][20]
    pos = list(set(cds_dict[ii].split(',')))
    pos = list(map(int, pos))
    pos.sort()
    pos = pos[pos.index(int(cds_dict[ii].split(',')[0])):(pos.index(int(cds_dict[ii].split(',')[1]))+1)]
    pos_list = getPosList(pos)
    ref_seq = seq[:]
    alt_seq = seq[:]
    exon_pos = []
    if annotation_dict[ii][2] == '+':
        for iii in range(len(position_info)):
            if position_info[iii][0:2] == 'ex':
                single_pos = int(subname[iii])
                exon_pos.append(single_pos)
                single_pos_index = pos_list.index(single_pos)
                ref_seq=ref_seq[:single_pos_index] + mnv_info[3].split(',')[iii] + ref_seq[single_pos_index+1:]
                alt_seq=alt_seq[:single_pos_index] + mnv_info[4].split(',')[iii] + alt_seq[single_pos_index+1:] 
    else: 
        pos_list = toReverse(pos_list)
        for iii in range(len(position_info)):
            if position_info[iii][0:2] == 'ex':
                single_pos = int(subname[iii])
                exon_pos.append(single_pos)
                single_pos_index = pos_list.index(single_pos)
                ref_seq=ref_seq[:single_pos_index] + codon_change[mnv_info[3].split(',')[iii]] + ref_seq[single_pos_index+1:]
                alt_seq=alt_seq[:single_pos_index] + codon_change[mnv_info[4].split(',')[iii]] + alt_seq[single_pos_index+1:]
        exon_pos = toReverse(exon_pos) 
    return [exon_pos,pos_list,ref_seq,alt_seq]   
# Obtain AA change type
def getCodonAAType(exon_pos,pos_list,ref_seq,alt_seq,codon2aa_dict):
    total = []
    codon = []
    aa_res = []
    aa_is_start=[]
    aa_types = []
    k = 0
    for i in range(len(exon_pos)):
        index_pos = pos_list.index(exon_pos[i])
        codon.append(ref_seq[index_pos] + str(index_pos+1) + alt_seq[index_pos])
        aa_pos = math.ceil((index_pos+1)/3)
        if i != 0:
            if aa_pos == math.ceil((pos_list.index(exon_pos[i-1])+1)/3):
                aa_res[k-1] = '@' + aa_res[k-1]
                aa_is_start[k-1] = '@' + aa_is_start[k-1]
            else:
                aa_ref = ref_seq[(aa_pos*3-3):(aa_pos*3)]
                aa_alt = alt_seq[(aa_pos*3-3):(aa_pos*3)]
                aa_res.append(codon2aa_dict[aa_ref] + str(aa_pos) + codon2aa_dict[aa_alt])
                aa_is_start.append('0')
                k = k+1
        else:
            aa_ref = ref_seq[(aa_pos*3-3):(aa_pos*3)]
            aa_alt = alt_seq[(aa_pos*3-3):(aa_pos*3)]
            aa_res.append(codon2aa_dict[aa_ref] + str(aa_pos) + codon2aa_dict[aa_alt])
            if pos_list.index(exon_pos[i])<=2:
                aa_is_start.append('1')
            else:
                aa_is_start.append('0')
            k = k+1
    for i in range(len(aa_res)):
        single_aa = aa_res[i].split('@')[-1]
        single_aa_is_start = aa_is_start[i].split('@')[-1]
        aa_types.append(getAAChangeType(single_aa[0],single_aa[-1],single_aa_is_start))
    aa_types = getUniqAndSort(aa_types)
    total.append(codon)
    total.append(aa_res)
    total.append(aa_types)
    return total     
# Adjust position_info
def fixPositionInfo(position_info):
    temp = []
    for x in position_info:
        temp.append(x.split(':')[0])
    return temp
# Obtain risk type：High，Moderate，Low，Modifier
def getRisk(types_list):
    temp = str(types_list)
    flag = 'Modifier'
    if 'synonymous_variant' in temp or 'start_gain' in temp or 'splice_region_variant' in temp:
        flag = 'Low' 
    if 'missense_variant' in temp:
        flag = 'Moderate'
    if 'stop_gain' in temp or 'stop_lost' in temp or 'start_lost' in temp or 'splice_donor_variant' in temp or 'splice_acceptor_variant' in temp: 
        return 'High'
    return flag
####################################################################

print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Run MNVAnnoGene.py')
########参数读入############################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load parameters')
mnv_file_url = args.infile
annotation_file_url = args.database
if not os.path.exists(annotation_file_url):
    print('no databases')
    exit()
chr_list = args.chr
upstream = int(args.upstream)
downstream = int(args.downstream)
splice = int(args.splice)
output = args.output
if '/' in mnv_file_url:
    mnv_result_url = '/'.join(mnv_file_url.split('/')[:-1])+'/'+output
else:
    mnv_result_url=output
chr_list=chr_list.split(',')
############################################################################

########Data loading, processing and output################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load data')
# Create dict
codon2aa_dict = createCodon()
types_dict = createTypes()
codon_change = createCodonChange()
# Load database
chr_anno_dict=readAllAnno(chr_list,annotation_file_url)
# Load MNV file
chr_mnv_dict=readAllMNV(chr_list,mnv_file_url)
# Open files
result=open(mnv_result_url,'w')
# Cycle for chr
for chrid in chr_list:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Process the chr' + chrid)
    # Extract single chr
    annotation_dict = readAnnotation(chr_anno_dict,chrid)
    if len(annotation_dict)==0:
        for i in chr_mnv_dict[chrid]:
            result.write(i.strip('\n')+'\t.\n')        
        continue    
    # Adjust database
    annotation_dict = fixAnnotation(upstream,downstream,splice,annotation_dict)
    # Create point_dict
    point_dict = createPoint(annotation_dict)
    # Create range_dict
    range_dict=createRange(annotation_dict)
    # Create cds_dict
    cds_dict = createCDS(annotation_dict)
    # Extract MNV for the chr
    mnv_list = readMNV(chr_mnv_dict,chrid)
    # Create mnv_dict
    mnv_dict = createMNV(mnv_list)
    # Annotate single point
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Annotating the information of each point of each MNV')
    mnv_dict = snvAnno(range_dict,mnv_dict,point_dict,annotation_dict)
    # Annotate MNV
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Annotating the information of each MNV and output')
    for i in mnv_list:
        # Extract each point information of MNV
        res = getTrans(i,mnv_dict)
        subname = res[0]
        subdict = res[1]
        trans = res[2]
        mnv_all_trans=[]
        # Cycle transcript
        for ii in trans:
            mnv_trans=[]
            position_info = []
            # Judge if MNV falls on exon
            flag_exon_utr5 = judgeExonUTR5(subname,ii,subdict)
            flag_utr5 = flag_exon_utr5[0]
            flag_exon = flag_exon_utr5[1]
            mnv_trans.append(annotation_dict[ii][3])
            mnv_trans.append(annotation_dict[ii][4])
            mnv_trans.append(ii)
            position_info = getInfo(subname,ii,subdict)
            mnv_trans.append(','.join(position_info))
            # Judge start_gain for MNV falling on 5'UTR
            if flag_utr5 !=0:
                res_utr5 = getSeqInfoUTR5(i,ii,annotation_dict,position_info,subname)
            # Judge start_gain for MNV falling on exon
            if flag_exon != 0 :
                res = getSeqInfo(i,ii, annotation_dict,cds_dict,position_info,subname)
                res = getCodonAAType(res[0],res[1],res[2],res[3],codon2aa_dict)
                if annotation_dict[ii][2] =='+':
                    mnv_trans.append(','.join(res[0]))
                    mnv_trans.append(','.join(res[1]))
                else:
                    mnv_trans.append(','.join(toReverse(res[0])))
                    mnv_trans.append(','.join(toReverse(res[1])))
            else:
                mnv_trans.append('.')
                mnv_trans.append('.')        
            position_info = fixPositionInfo(position_info)
            mm=getUniqAndSort(position_info)
            types = {x:types_dict[x] for x in mm if x in types_dict}
            # Add types of 5'UTR and exon
            if flag_utr5 != 0 :
                types_list = list(types.values()) + [res_utr5]
            else:
                types_list = list(types.values())
            if flag_exon != 0 :
                if annotation_dict[ii][2] =='+':
                    types_list = types_list + res[2]
                else:
                    types_list = types_list + toReverse(res[2])
            mnv_trans.append(','.join(types_list))
            # Add risk types
            mnv_trans.append(getRisk(types_list))
            # Splice
            mnv_all_trans.append(' '.join(mnv_trans))
        # Splice again
        temp = i[0:6]
        if mnv_all_trans:
            temp.append('|'.join(mnv_all_trans))
        else:
            temp.append('.')
        result.write('\t'.join(temp)+'\n')   
# Close files
result.close()
############################################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'End MNVAnnoGene.py')






