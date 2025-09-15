## 比较了不同人群的MNV差异。【有哪些维度进行比较】
## 不同人群的MNV重叠情况 ===================
## 统同不同人群的韦恩图
    python 01getMNVVenn.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv ./venn/venn_1000G.txt
## 因为人群间unique的占比很高，所以为了找到原因做了两项工作
## 查看SNV的人群韦恩图
    # 获得不同人群的SNV数据
        python 02get1000GSNV.py /data/jinww/04reference/publicDB/1000G/hg38/ snv/ 
    # 统计不同人群的SNV的overlap情况
        python 03getSNVVenn.py获得uniqu的SNV的概率，发现之所以韦恩图出现上面的原因，是因为SNP的分布导致
    # 查看SNV的不同人群的密度,5个人群都鉴定
        python 04getDes.py snv/AFR.snv snv/AFR_snv_density 
## 将不同人群共有的SNV的MNV提取出来，然后根据这个绘制真正的不同人群的韦恩图
    python 05getMNVfromSameSNV.py  /data/jinww/mnv/analyse3/02pop/01overlap/ /data/jinww/mnv/analyse3/02data/01origin/hg38mnv
    # 获得common_snv 和venn_1000G_common
## 获得不同人群的MNV密度 (这里不是共有的，是所有的)
    python 06getPopMNV.py /data/jinww/04reference/publicDB/1000G/EAS/all.EAS.adjust.sorted.mnv EAS EAS.mnv
    python 04getDes.py mnv/EAS.mnv mnv/EAS_mnv_density

