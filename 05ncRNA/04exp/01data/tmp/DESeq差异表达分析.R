library(DESeq2)
d<-read.table('表达矩阵_SNV1_NC.txt',header=TRUE)
coldata<-read.table('分组文件_SNV1_NC.txt',header=TRUE)
dds <- DESeqDataSetFromMatrix(countData = d, colData = coldata, design= ~ condition) 
dds <- DESeq(dds)
res <- results(dds)
table(res$padj <0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "差异表达数据_SNV1_NC.csv")

