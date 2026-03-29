import sys

### MHC chr6:28477797-33448354 (hg19)                chr6:28510120-33480577 (hg38)
fo = open(sys.argv[2],"w")
with open(sys.argv[1],"r") as fi:
    fo.write(fi.readline())
    for line in fi:
        data = line.strip().split()
        if sys.argv[3] == "hg19":
            if data[0] == "6" and int(data[2]) in range(28477797,33448355):
                continue
            else:
                fo.write(line)
        elif sys.argv[3] == "hg38":
            if data[0] == "6" and int(data[2]) in range(28510120,33480578):
                continue
            else:
                fo.write(line)
fi.close()
fo.close()
