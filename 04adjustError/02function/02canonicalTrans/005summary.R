library(data.table)
library(stringr)
file=c('gnomAD','fiveData','total')
for (filename in file) {
  print(filename)
  filename=paste0('merge_',filename,'.txt')
  df<-read.table(filename,header = T,sep='\t',stringsAsFactors = F)
  df<-df[c(3,7,9,15,19,23,27)]
  df[df$snv3AAC=='.','snv3AAC']<-'K/AA'
  df$MNVAAC<-str_split_fixed(df$MNVAAC,'[/]',n=Inf)[,2]
  df$snv1AAC<-str_split_fixed(df$snv1AAC,'[/]',n=Inf)[,2]
  df$snv2AAC<-str_split_fixed(df$snv2AAC,'[/]',n=Inf)[,2]
  df$snv3AAC<-str_split_fixed(df$snv3AAC,'[/]',n=Inf)[,2]
  judgeAA<-function(x){
    if(x[1]!=x[2] & x[1]!=x[3] & x[1]!=x[4]){
      return(1)
    }else{
      return(0)
    }
  }
  df$flag<-apply(df[,4:7],1,judgeAA)
  df_complete<-df[df$flag==1,]
  print(paste('mnvid:',length(unique(df_complete$MNVID))))
  print(paste('gene:',length(unique(df_complete$gene))))
  print(paste('transcript:',length(unique(df_complete$trans))))
}
