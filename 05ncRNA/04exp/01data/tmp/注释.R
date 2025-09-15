#setwd("D:/Desktop/202009小鼠转录组/9天/lncRNA/differential_expression_analysis_SMKI/htseq")
gtf=read.table("Homo_gtf_information.txt",sep="\t",head=F)
colnames(gtf)=c("Row.names"," gene_name","chr","start","end","forword","kind")
head(gtf)
gtf1=gtf[,c("Row.names"," gene_name","chr","start","end","forword","kind")]
head(gtf1)

data=read.csv("差异表达数据_SNV2_NC.csv",sep=",",head=T)
head(data)

df=merge(gtf1,data,by="Row.names",all.y=T)
head(df)

write.table(df,file = "differential_expression_gene_SNV2_NC.txt",sep = "\t",row.names=FALSE,quote=F)
