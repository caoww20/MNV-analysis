library(data.table)
pop<-c('AFR','AMR','EAS','EUR','SAS')
# Read all files
combinations <- combn(pop, 2)  # Generate all pairwise combinations
combinations<-apply(combinations, 2, paste0,collapse='_')
df<-fread('./03fst/filter/AFR_AMR.weir.fst',header = T,data.table = F)
df<-df[1:2]
subpop_list=list()
for (i in combinations) {
  tmp<-fread(paste0('./03fst/filter/',i,'.weir.fst'),header = T,data.table = F)
  subpop_list[[i]]<-tmp$WEIR_AND_COCKERHAM_FST
}
res<-cbind(df,do.call(cbind,subpop_list))
# Count statistics
judgeSig<-function(x){
  x<-x[!is.na(x)]
  f1=0
  f2=0
  flag<-sum(x>0.15)
  if (flag>=1) {
    f1=1
  }
  if (flag==4){
    f2=1
  }
  return(c(f1,f2))
}
judgeSig2<-function(x){
  x<-x[!is.na(x)]
  f1=0
  flag<-sum(x>0.15)
  if (flag==4){
    f1=1
  }
  return(f1)
}
getNum<-function(df,key_pop){
  flag<-apply(df[,grepl(key_pop,colnames(df))], 1, judgeSig)
  return(c(key_pop,rowSums(flag)))
}
getDf<-function(df,key_pop){
  newdf<-cbind(df[1:2],df[,grepl(key_pop,colnames(df))])
  flag<-apply(newdf[3:6], 1, judgeSig2)
  return(newdf[flag==1,1:2])
}
res_list=list()
for (keypop in pop) {
  res_list[[keypop]]=getNum(res,keypop)
  write.table(getDf(res,keypop),paste0('03fst/significant/sig',keypop,'.txt'),row.names = F,sep = '\t',quote = F)
}
res2<-as.data.frame(do.call(rbind,res_list))
colnames(res2)<-c('pop','total','sig')
res2<-reshape2::melt(res2,id=c('pop'))
res2$value<-as.numeric(res2$value)
write.table(res2,'03fst/significant/sig.txt',row.names = F,sep = '\t',quote = F)


