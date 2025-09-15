# 准备工作
    Homo_sapiens.GRCh38.108.chr.gtf Homo_sapiens.GRCh38.cdna.all.fa chr_fix.fa # 仅包含常见染色体和MT
    cut -f1 -d. Homo_sapiens.GRCh38.cdna.all.fa > cdna.fa
    awk -F '\t' 'BEGIN{OFS="\t"} {print $2,$1}' /data/jinww/mnv/library/human/01gtf/ensembl/hg38/transcript.list >trans.map
    awk -F '\t' 'BEGIN{OFS="\t"} {print $1,$4,$5}' /data/jinww/mnv/MNVAnno/database/human/en_GRCh38.108_gene.anno.txt|sort -u >symbol.map
    # 创建索引文件夹
    mkdir bowtie salmon_mRNA
    # salmon
        grep "^>" chr_fix.fa | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt
        cat Homo_sapiens.GRCh38.cdna.all.fa chr_fix.fa > gentrome_mRNA.fa
        salmon index -t gentrome_mRNA.fa -d decoys.txt -p 10 -i salmon_mRNA/
    # salmon (版本问题，所以不能有参数-d)
        salmon index -t cdna.fa -p 40 -i salmon_mRNA/
    # bowtie
        bowtie-build chr_fix.fa bowtie/genome
# mRNA
## 清洗数据
    cat sample.list|while read id;do (fastp -i ${id}_1.fastq -I ${id}_2.fastq -o ${id}_1.clean.fq -O ${id}_2.clean.fq -w 5 -q 20 -h ${id}_clean.html &);done
    # -w 是线程 -q 是质量 fastp自动去除接头
## 基因定量 （salmon）
    salmon quant -i ../../../index/salmon_mRNA/ -l IU -g ../../../index/Homo_sapiens.GRCh38.108.chr.gtf  -1 02T0001_1.clean.fq -2 02T0001_2.clean.fq -p 20 -o ../03salmon/02T0001_quant &
    cat ../01data/sample.list|while read id;do (salmon quant -i ../../../index/salmon_mRNA/ -l IU -g ../../../index/Homo_sapiens.GRCh38.108.chr.gtf -1 ${id}_1.clean.fq -2 ${id}_2.clean.fq -p 20 -o ../03salmon/${id}_quant &);done
    # -l A 表示自动识别测序类型 # -g 基因定量
    # 由于版本问题A无法使用，所以改用-l IU # IU表示无连特异性
## 合并结果   
    Rscript 001merge.R # 这个是在代码里运行结果

# miRNA定量
## 已知的miRNA下载
    # miRBase
    mature.fa.gz hairpin.fa.gz gff3
    ## 提取mature
    extract_miRNAs.pl mature.fa hsa > hsa_mature.fa
    ## 提取hairpin
    extract_miRNAs.pl hairpin.fa hsa > hsa_hairpin.fa
## 清洗数据
    cat sample.list|while read id;do (fastp -i ${id}.fastq -o ${id}.clean.fq -w 5 -q 20 -h ${id}_clean.html &);done
## 使用mirdeep2进行定量（太慢，使用自用代码）
    ### 比对
        cat sample.list|while read id;do (mkdir m_${id});done
        cat sample.list|while read id;do (mapper.pl ${id}.clean.fq -e -h -j -l 18 -o 20 -m -p ~/04reference/zebrafish/bowtie/genome -s m_${id}/${id}_reads.fa -t m_${id}/${id}_reads.arf &);done
    ### 定量
        cat sample.list|while read id;do (quantifier.pl -p ~/05public/miRBase/dme_hairpin.fa -m ~/05public/miRBase/dme_mature.fa -r m_${id}/${id}_reads.fa -t hsa -y now 2> Quantifier_log &);done 
        # 如果有新的则需要更改

    bowtie-align --wrapper basic-0 -p 20 -f -n 0 -e 80 -l 18 -a -m 5 --best --strata ../../../index/bowtie/genome --al dir_mapper_seq_01S0015.clean.fq_2124579877_04_12_2023_t_22_17_50/01S0015.clean.fq_mapped --un dir_mapper_seq_01S0015.clean.fq_2124579877_04_12_2023_t_22_17_50/01S0015.clean.fq_not_mapped dir_mapper_seq_01S0015.clean.fq_2124579877_04_12_2023_t_22_17_50/reads_nr.fa dir_mapper_seq_01S0015.clean.fq_2124579877_04_12_2023_t_22_17_50/mappings.bwt
## 由于mirdeep2进行定量特别慢，所以使用自己的代码
    # python 001quantifier.py mature_seq.list SRR.fa SRR.quant
    cat 01data/sample.list |while read id;do (python 001quantifier.py seq.list 02clean/${id}.clean.fq 03quant/${id}.quant &);done
## 合并结果
    Rscript 002merge.R # 这个是在代码里运行结果

## 绘制图片
    miRNA的堆叠图
    KEGG分析
    韦恩图，唯一差异的gene
    具体的demo例子