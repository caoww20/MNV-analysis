df<-read.table('mnv.genotype',header = F,sep = '\t',stringsAsFactors = F)
df$id<-apply(df[c(2,4,5)],1,paste,collapse='-')
df<-df[1044:nrow(df),]
tmp<-df[8:32]
tmp[tmp >=1]<-1
tmp$sum<-apply(tmp,1,sum)
res<-cbind(tmp$sum,df$id)
colnames(res)<-c('num','flag')
write.table(res,'Qing_SN.txt',sep = '\t',quote = F,row.names = F)

tmp<-df[8:32]
tmp$sum<-apply(tmp,1,sum)
res<-cbind(tmp$sum,df$id)
colnames(res)<-c('AC','flag')
write.table(res,'Qing_AC.txt',sep = '\t',quote = F,row.names = F)


# 验证我提取的样本是否正确
ref<-read.table('mnv.txt',sep = '\t',stringsAsFactors = F)
res2<-cbind(ref[1044:nrow(ref),9],tmp$sum)
res2<-as.data.frame(res2)
unique(res2$V2-res2$V1)
# 为0 是正确的