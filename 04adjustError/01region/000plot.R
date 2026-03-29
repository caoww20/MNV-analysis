library(ggplot2)
##主题
{
  mythem<-theme(panel.grid=element_blank(),
                plot.margin=unit(rep(2,4),'lines'),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                title = element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.text.x=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black')
  )
  mythem2<-theme(panel.grid=element_blank(),
                plot.margin=unit(rep(2,4),'lines'),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                title = element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.text.x=element_text(vjust=0.5,size=16,angle = 45,face = "bold", color='black'),
                axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black')
  )
}

# 各个数据集的区域分布情况
{
  fun1<-function(f){
    df<-read.table(f,header = F,stringsAsFactors = F,sep = '\t')
    colnames(df)<-c('type','num')
    return(df)
  }
  df<-cbind(fun1('01num/GnomAD.txt'),fun1('01num/1000G.txt'),fun1('01num/GTEx.txt'),fun1('01num/UKB50w.txt'),fun1('01num/TCGA.txt'),fun1('01num/UKB20w.txt'))
  df<-df[c(1:2,4,6,8,10,12)]
  colnames(df)<-c('type','GnomAD','1000G','GTEx','UKB50w','TCGA','UKB20w')  
  df$id<-1:nrow(df)
  res<-reshape2::melt(df,id=c('type','id'))  
  
  p<-ggplot(res,aes(x=id,y=log10(value),fill=variable))+
    geom_bar(stat = 'identity',position = "dodge",show.legend = T)+
    # geom_bar(stat = 'identity',show.legend = T)+
    labs(title='',x='Annotation Type',y='Number of MNV (log10)')+
    scale_x_continuous(breaks=c(1:nrow(df)),labels=df$type)+
    theme_bw()+mythem+
    scale_fill_manual(values = c("GnomAD"="black","1000G" = "#1A325F", "GTEx" = "#2A5888", "UKB50w" = "#4396B1", "TCGA"="#89CEED","UKB20w"="#C4F5FC"))
  p
  
  ggsave(
    filename = './graph/bar.pdf',
    plot = p,width = 8,height = 3,
    units = 'in'
  )
  
}

# 各个数据集的区域密度分布情况
# density
{
  fun1<-function(f){
    df<-read.table(f,header = F,stringsAsFactors = F,sep = '\t')
    df$V2<-df$V2*(1e+06)/mnvRegion$V2
    colnames(df)<-c('type','freq')
    df$id<-1:nrow(df)
    return(df)
  }
  pgraph<-function(df,f,col1,a,b){
    p<-ggplot(df, aes(x = id,y=freq)) +
      geom_point(size=2,color=col1)+
      geom_line(lwd=1,color=col1) +
      labs(title='',x='',y='Density of MNVs')+
      scale_x_continuous(breaks=c(1:nrow(df)),labels=df$type)+
      # ylim(0,1200)+
      scale_y_continuous(breaks=a, limits = b)+  # ylim 和scale_y_continuous不能同时用
      theme_bw()+mythem
    ggsave(
      filename = paste0('./graph/',f,'.pdf'),
      plot = p,width = 4,height = 3,
      units = 'in'
    )
  }
  mnvRegion<-read.table('01num/region.txt',header = F,stringsAsFactors = F,sep = '\t')
  G1000<-fun1('01num/1000G.txt')
  pgraph(G1000,'1000G','#1A325F',seq(0, 1500, length.out = 4),c(0,1500))

  GTEx<-fun1('01num/GTEx.txt')
  pgraph(GTEx,'GTEx','#2A5888',seq(0, 600, length.out = 4),c(0,600))
  
  UKB50w<-fun1('01num/UKB50w.txt') 
  pgraph(UKB50w,'UKB50w','#4396B1',seq(0, 450, length.out = 4),c(0,450))
  
  TCGA<-fun1('01num/TCGA.txt')
  pgraph(TCGA,'TCGA','#89CEED',seq(0, 90, length.out = 4),c(0,90))
  
  UKB20w<-fun1('01num/UKB20w.txt')
  pgraph(UKB20w,'UKB20w','#C4F5FC',seq(0, 900, length.out = 4),c(0,900))
  
  gnomAD<-fun1('01num/gnomAD_WGS.txt')
  pgraph(gnomAD,'gnomAD','black',seq(0, 2400, length.out = 4),c(0,2400))

}
