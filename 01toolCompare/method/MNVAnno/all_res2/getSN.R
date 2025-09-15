df<-read.table('mnv.genotype',header = F,sep = '\t',stringsAsFactors = F)
df$id<-apply(df[c(2,4,5)],1,paste,collapse='-')
tmp<-df[8:32]
tmp[tmp >=1]<-1
tmp$sum<-apply(tmp,1,sum)
res<-cbind(tmp$sum,df$id)
colnames(res)<-c('num','flag')
write.table(res,'MNVAnno_SN.txt',sep = '\t',quote = F,row.names = F)
