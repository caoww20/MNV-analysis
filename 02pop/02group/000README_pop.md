## MNV对于人群分层的影响 (使用的是老数据，因为新数据增加了MHC、black和X，但不影响这里)
## 拿共有的SNV的MNV作为研究对象，看其在人群分层和进化等方面的研究。【聚类效果查看】
## 初始数据获得
    ## 先从luohh的文件夹下拼接1000G的mnv和mnv_multi
    ## 获得样本信息
    python 001getPopSample.py /data/jinww/04reference/publicDB/1000G/igsr-1000_genomes_30x_on_grch38.tsv  /data/jinww/mnv/analyse3/02pop/02group/01origin/sample
## 共有SNV的MNV捕获
    # 从共有的SNV中获得MNV 结果为MNVFromCommonSNV.txt
        python 001getMNVFromCommonSNV.py /data/jinww/mnv/analyse3/02pop/01overlap/venn/venn_1000G_common.txt /data/jinww/mnv/analyse3/02data/02datasets/1000G.txt 01origin/mnv.txt 02adjust/MNVFromCommonSNV.txt &
    # 获得基因型
        python /data/jinww/mnv/MNVAnno/tools/MNVGetGenotype.py -m MNVFromCommonSNV.txt -v /data/jinww/04reference/publicDB/1000G/hg38/IGSR.chr1.vcf -g 1 &
    # 将其转化为vcf格式，并且对位置重叠的MNV进行调整即+n
        # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
        # 1 85  mnv1    A   T   .   PASS    .   GT
        python 002toVCF.py 01origin/VCF.head 02adjust/MNVFromCommonSNV.genotype 02adjust/MNVFromCommonSNV.vcf &
## 共有MNV捕获
    python 003getCommonMNV.py /data/jinww/mnv/analyse3/02pop/01overlap/venn/venn_1000G_common.txt 02adjust/MNVFromCommonSNV.vcf 02adjust/commonMNV.vcf &
## 捕获频率 # 004getAF.py比较慢，改成版本2了
    python 004getAF2.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv 02adjust/MNVFromCommonSNV.vcf 02adjust/MNVFromCommonSNV.freq
    python 004getAF2.py /data/jinww/mnv/analyse3/02data/01origin/hg38mnv 02adjust/commonMNV.vcf 02adjust/commonMNV.freq
## 进行质控，剔除高缺失率（--geno 0.05）和极低等位基因频率( --maf 0.01 )的SNP 1000G基因组的数据，使用0.001，暂时不考虑罕见变异
    plink2 --vcf MNVFromCommonSNV.vcf --geno 0.05 --maf 0.01 --recode vcf --out MNVFromCommonSNVFilter &
    plink2 --vcf commonMNV.vcf --geno 0.05 --maf 0.01 --recode vcf --out commonMNVFilter &
    ## 1000G.vcf软连接到这里 # 已经做了质控
## 计算PCA, 需要到vcf的文件处 【做过3种机制的PCA，但差不多，所有就不用了】
    bash ../006plink.sh MNVFromCommonSNVFilter &
    bash ../006plink.sh commonMNVFilter &
    bash ../006plink.sh 1000G & 
## 计算fst
    python 007getFSTSH.py MNVFromCommonSNVFilter.vcf AFR,AMR,EAS,EUR,SAS a.sh
    python 007getFSTSH.py 1000G.vcf AFR,AMR,EAS,EUR,SAS a.sh
    bash a.sh
    # 单点的结果
        CHROM：染色体
        POS：SNP位置
        FST：单个位点Fst
    # 区间的结果
        CHROM: 表示染色体编号，这里为1。
        BIN_START: 基因组区间的起始位置。
        BIN_END: 基因组区间的结束位置。
        N_VARIANTS: 该基因组区间内的变异位点数量。
        WEIGHTED_FST: 加权 Fst 值，是一种考虑变异位点频率差异的群体遗传分化指数。
        MEAN_FST: 平均 Fst 值，是一种简单平均计算的群体遗传分化指数
        加权 Fst 考虑了变异位点频率差异，对于频率差异较大的位点给予更高的权重，更准确地反映了种群间的基因差异。这在研究不同种群之间的遗传分化程度时可能更有意义。
        平均 Fst 是简单平均计算的群体遗传分化指数，对于数据中的所有位点一视同仁地考虑，可能更容易计算和解释。在样本多样性较低、变异位点频率差异较小的情况下，平均 Fst 也可以提供一定的研究价值。
        综上所述，根据研究目的和数据特征，可以灵活选择使用加权 Fst 或平均 Fst，或者结合两者综合考虑。
    # FST取值范围的意义
        Fst取值范围[0,1]，
        Fst=0时，表明不同群体遗传结构完全一致，群体间没有分化，
        Fst=1时，表明等位基因在不同的群体中被固定，完全分化。
        Wright提出，在实际研究中，
        如果Fst取值为0～0.05群体间遗传分化很小，可忽略不计；
        Fst取值为0.05～0.15，群体间存在中等程度的遗传分化；
        Fst取值为0.15～0.25，群体间遗传分化较大；
        Fst取值为>0.25时，群体间有很大的遗传分化。
    # 其他备注信息
        # 群体间遗传分化指数fst的结果是负数是什么意思
        # 群体间遗传分化指数 Fst 是用来衡量不同种群之间基因差异的指标。如果 Fst 的结果是负数，通常情况下这意味着不同种群间的基因差异比总体内个体间的基因差异要小。这种情况可能是由于种群间的基因交流或者混合造成的。需要注意的是，Fst 的结果还受到数据质量、采样策略等多种因素的影响，因此需要谨慎地解释结果
        # 当群体间遗传分化指数 Fst 的结果为 "NaN" 时，这通常表示计算过程中出现了非数值或未定义的情况，导致无法得到合理的结果。这可能是由于样本数量不足、基因数据缺失或其他技术问题导致的。需要对数据质量和分析方法进行仔细检查，以确定如何进一步处理数据以获得可靠的结果。
        # 过滤掉负的点和区域
    ## 合并结果
        python 008mergeFST.py ./03fst/ AFR,AMR,EAS,EUR,SAS 0 five_single_merge &  # single只会得到一个结果
        python 008mergeFST.py ./03fst/ AFR,AMR,EAS,EUR,SAS 1 five_region_merge &  # region会得到两个结果，一个是mean一个是weight
    ## 这里我也做了没有filter（即使用MNVFromCommonSNV.vcf文件）的结果，发现其实差不多
## fst可以计算出在人群中有分化的位点 这些点意味着进化
## 【舍弃】 提取分化指数0.05 0.15 0.2的的数据提取vcf文件并作pca分析，但是结果和之前的pca差不多，所以舍弃
## 【舍弃】 上述的结果也进行富集分析查看，没有显著的通路（也肯能是我不怎么接触），主要是基因的数量也特别多，所以舍弃
    
## 计算不同人群的FST≥0.15的MNV的数量（以及任一都大于0.15的数量）
    Rscript 011getSigNum.R # 得到sig.txt和各个人群的sig数据
## 【舍弃】 AFR为例，看一下任一都大于0.15的数量对应的基因位于哪个上面 # 没有得到有意思的结果，虽然富集到心肌病中
## 使用all.weir.fst获取SNV和MNV的FST≥0.15的数据，并提取基因内至少包含10个及以上FST显著MNV，进行富集分析
    awk '$3 >= 0.15 {print $1"\t"$2"\t"".\t""A\tT\t."}' 03fst/snv/all.weir.fst|grep -v CHROM|sed 's/chr//g' > 03fst/significant/snv_sig.txt &
    bash /data/jinww/mnv/MNVAnno/MNVAnno.sh -i snv_sig.txt -t gene -D en_GRCh38.108_gene.anno.txt -o snv_sig.gene
    awk '$3 >= 0.15 {print}' 03fst/filter/all_5.weir.fst > 03fst/significant/mnv_sig.txt &
    python 012getInfo2.py 03fst/significant/mnv_sig.txt 02adjust/MNVFromCommonSNVFilter.vcf /data/jinww/mnv/analyse2/02data/02annotation/gene.txt 03fst/significant/mnv_sig.gene &
    Rscript 013getGene.R # 获得基因对应的FST大于0.15的数量，取包含10个以上FST显著的做富集分析
    # 使用david工具进行富集分析
    # 因为SNV的数量已经多到没法衡量，所以只考虑MNV，最终结果是MNV_GO.txt


## 获得cardiomyopathy的数据
    python 014getcardiomyopathyMNV.py dilated_cardiomyopathy.gwas /data/jinww/mnv/analyse3/02data/02datasets/1000G.txt overlap.MNV


