file_list <- list.files(full.names = TRUE)
trans_tpm=list()
trans_count=list()
gene_tpm=list()
gene_count=list()
# 读入文件
k=0
for (i in file_list) {
  id=gsub('./','',i)
  id=gsub('_quant','',id)
  if (k==0) {
    trans_url=paste0(i,'/quant.sf')
    tmp<-read.table(trans_url,header = T,stringsAsFactors = F,sep = '\t')
    trans_tpm[[id]]<-cbind(tmp[1:3],round(tmp[c(4)],2))
    trans_count[[id]]<-cbind(tmp[1:3],round(tmp[c(5)],2))
    gene_url=paste0(i,'/quant.genes.sf')
    tmp<-read.table(gene_url,header = T,stringsAsFactors = F,sep = '\t')
    gene_tpm[[id]]<-cbind(tmp[1:3],round(tmp[c(4)],2))
    gene_count[[id]]<-cbind(tmp[1:3],round(tmp[c(5)],2))
    k=k+1
  }else{
    trans_url=paste0(i,'/quant.sf')
    tmp<-read.table(trans_url,header = T,stringsAsFactors = F,sep = '\t')
    trans_tpm[[id]]<-round(tmp[c(4)],2)
    trans_count[[id]]<-round(tmp[c(5)],2)
    
    gene_url=paste0(i,'/quant.genes.sf')
    tmp<-read.table(gene_url,header = T,stringsAsFactors = F,sep = '\t')
    gene_tpm[[id]]<-round(tmp[c(4)],2)
    gene_count[[id]]<-round(tmp[c(5)],2)
  }
}
sample=c()
for (i in file_list) {
  id=gsub('./','',i)
  id=gsub('_quant','',id)
  sample=c(sample,id)
}

# 合并文
trans_tpm=do.call(cbind,trans_tpm)
colnames(trans_tpm)<-c('trans','length','effect_length',sample)
trans_count=do.call(cbind,trans_count)
colnames(trans_count)<-c('trans','length','effect_length',sample)

gene_tpm=do.call(cbind,gene_tpm)
colnames(gene_tpm)<-c('gene','length','effect_length',sample)
gene_count=do.call(cbind,gene_count)
colnames(gene_count)<-c('gene','length','effect_length',sample)

# 仅保留TPM大于1的行
flag<-apply(trans_tpm[4:ncol(trans_tpm)],1,max)
trans_tpm<-trans_tpm[which(flag>1),]
trans_count<-trans_count[which(flag>1),]
                     
flag<-apply(gene_tpm[4:ncol(gene_tpm)],1,max)
gene_tpm<-gene_tpm[which(flag>1),]
gene_count<-gene_count[which(flag>1),]
# 将结果输出
write.table(trans_tpm,'trans_tpm.txt',quote = F,sep = '\t',row.names = F)
write.table(trans_count,'trans_count.txt',quote = F,sep = '\t',row.names = F)
write.table(gene_tpm,'gene_tpm.txt',quote = F,sep = '\t',row.names = F)
write.table(gene_count,'gene_count.txt',quote = F,sep = '\t',row.names = F)
