import sys

mnvInfoDict = {}
with open("/home/luohh/MNVAnalysis/01.GWAS/01.Height/01.Data/02.MNVData/UKB50w.filterMulti.2joint.autosome.mnv","r") as fi:
    for line in fi:
        data = line.strip().split()
        mnvInfoDict[data[2]] = data
fi.close()

sigMNVsRsidDict = {}
sigMNVIDDict = {}
sigMNVsFormatDict = {}
with open(sys.argv[1],"r") as fi:
    for line in fi:
        data = line.strip().split()
        rsidList = mnvInfoDict[data[1]][5].split(",")
        if "." in rsidList:
            continue
        else:
            sigMNVIDDict[data[1]] = 1
            sigMNVsFormatDict[data[1]] = "\t".join(mnvInfoDict[data[1]])
            for rsid in rsidList:
                sigMNVsRsidDict[rsid] = 1
fi.close()

fo = open(sys.argv[2],"w")
for mnvid in sigMNVIDDict:
    fo.write(mnvid+"\n")
fo.close()

fo = open(sys.argv[3],"w")
for rsid in sigMNVsRsidDict:
    fo.write(rsid+"\n")
fo.close()

fo = open(sys.argv[4],"w")
for key,value in sigMNVsFormatDict.items():
    fo.write(value+"\n")
fo.close()