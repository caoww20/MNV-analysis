file_list <- list.files(pattern = 'quant',full.names = TRUE)
res=list()
sample=c()
# 读入文件
k=0
for (i in file_list) {
  id=gsub('./','',i)
  id=gsub('.quant','',id)
  sample=c(sample,id)
  if (k==0) {
    tmp<-read.table(i,header = F,stringsAsFactors = F,sep = '\t')
    res[[id]]<-cbind(tmp[1:2],as.numeric(tmp$V3))
    k=k+1
  }else{
    tmp<-read.table(i,header = F,stringsAsFactors = F,sep = '\t')
    res[[id]]<-as.numeric(tmp$V3)
  }
}

# 合并文
res=do.call(cbind,res)
colnames(res)<-c('miRNA','seq',sample)
# 将结果输出
write.table(res,'miRNA_count.txt',quote = F,sep = '\t',row.names = F)

# 转化为RPM
res<-read.table('miRNA_count.txt',header = T,stringsAsFactors = F,sep = '\t')
rownames(res)<-res$miRNA
res<-res[3:ncol(res)]
read<-read.table('sample.read',header = F,sep = '\t',stringsAsFactors = F)
res<-rbind(res,read$V2)
df <- apply(res, 2, function(x) x*1000000 / x[length(x)])
df<-round(df,0)
write.table(df,'miRNA_RPM.txt',quote = F,sep = '\t',row.names = F)
