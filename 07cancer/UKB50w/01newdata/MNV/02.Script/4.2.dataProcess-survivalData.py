cancerList = ["01.Breast","02.Prostate","03.Melanoma","04.Colon","05.Lung","06.Cervical","07.Uterine","08.Bladder","09.Rectal","10.Ovary","11.Kidney"]

for cancer in cancerList:
    fo = open("/home/luohh/MNVAnalysis/01.GWAS/02.Cancer/01.Data/02.PhenoData/"+cancer+"/clinic.geno.012.tsv","w")

    clinicDict = {}
    with open("/home/luohh/MNVAnalysis/01.GWAS/02.Cancer/01.Data/02.PhenoData/"+cancer+"/clinic.tsv","r") as fi:
        header = fi.readline().split()
        fo.write("\t".join(header+["MNVID","genotype"])+"\n")
        for line in fi:
            data = line.strip().split()
            clinicDict[data[0]] = data
    fi.close()

    with open("/home/luohh/MNVAnalysis/01.GWAS/02.Cancer/03.Results/"+cancer+"/01.GenoData/01.MNVs/sigMNVs.vcf","r") as fi:
        for line in fi:
            if line[:2] == "##":
                continue
            elif line[:4] == "#CHR":
                sampleIDList = line.strip().split()
            else:
                data = line.strip().split()
                for key,value in clinicDict.items():
                    i = sampleIDList.index(key)
                    if data[i] == "0/0":
                        fo.write("\t".join(clinicDict[sampleIDList[i]]+[data[2],"0"])+"\n")
                    elif data[i] == "0/1":
                        fo.write("\t".join(clinicDict[sampleIDList[i]]+[data[2],"1"])+"\n")
                    elif data[i] == "1/0":
                        fo.write("\t".join(clinicDict[sampleIDList[i]]+[data[2],"1"])+"\n")
                    elif data[i] == "1/1":
                        fo.write("\t".join(clinicDict[sampleIDList[i]]+[data[2],"2"])+"\n")
                    else:
                        fo.write("\t".join(clinicDict[sampleIDList[i]]+[data[2],"NA"])+"\n")
    fi.close()
    