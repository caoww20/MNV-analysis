df<-read.table('03fst/significant/mnv_sig.gene.oneline.txt',header = F,stringsAsFactors = F,sep = "\t")
df$gene<-str_split_fixed(df$V7,' ',n=Inf)[,2]
df<-df[df$gene!='.',c(3,8)]
df<-df[!duplicated(df),]
gene<-as.data.frame(table(df$gene))
gene<-gene[order(-gene$Freq),]
# gene<-gene[gene$Freq>=10,]
write.table(gene,'03fst/significant/mnv_sig.gene.list',row.names = F,col.names = F)


df<-read.table('03fst/significant/snv_sig.gene.oneline.txt',header = F,stringsAsFactors = F,sep = "\t")
df$gene<-str_split_fixed(df$V7,' ',n=Inf)[,2]
df$id<-paste0(df$V1,'_',df$V2)
df<-df[df$gene!='.',c(9,8)]
df<-df[!duplicated(df),]
gene<-as.data.frame(table(df$gene))
# 提取前10%的基因进行富集分析
gene<-gene[order(-gene$Freq),]
# gene<-gene[1:as.integer(0.01*nrow(gene)),]
write.table(gene,'03fst/significant/snv_sig.gene.list',row.names = F,col.names = F)
