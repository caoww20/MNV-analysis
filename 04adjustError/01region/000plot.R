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

# 获得保守性得分图
{
  filterExtremum<-function(x){
    # x <- c(23, 17, 14, 29, 37, 21, 16, 18, 26, 31)
    iqr <- IQR(x)
    upper_bound <- median(x) + 1.5 * iqr
    lower_bound <- median(x) - 1.5 * iqr
    x_filtered <- x[x >= lower_bound & x <= upper_bound]
    return(x_filtered)
  }
  
  df<-read.table('02Cos_GTF/03pos_cos/res.cos',header = F,stringsAsFactors = F,sep = '\t')
  colnames(df)<-c('type','chr','pos','value')
  df<-na.omit(df)
  df$type<-factor(df$type)
  
  library(dplyr)
  # 假设你的数据框为df，包含"group"和"value"两列
  result_mean <- df %>%
    group_by(type) %>%
    summarize(avg_value = mean(value))
  write.table(result_mean,'graph/cos_average.txt',row.names = F,sep = '\t',quote = F)
  # result_median的结果不好，不使用
  result_median <- df %>%
    group_by(type) %>%
    summarize(avg_value = median(value))
  
  
  # 假设你的向量为vec
  res_list=list()
  for (i in c('intergenic','intro','atac','enhancer','lnc','piR','circ','miRNA','utr3','utr5','miRBS','TFBS','cds','splice','CE'  )) {
    categories <- cut(df[df$type==i,'value'], breaks = seq(0, 1, 0.1), include.lowest = TRUE, labels = FALSE)
    res_list[[i]]=table(categories)
  }
  res<-do.call(rbind,res_list)
  write.table(t(res),'graph/cos_group.txt',sep = '\t',quote = F)
  res_divided <- apply(res,1, function(x) x/sum(x))
  write.table(res_divided,'graph/cos_group_freq.txt',sep = '\t',quote = F)
  # 按照第10个模块的大小排序
  res_divided<-t(res_divided)
  res_divided<-res_divided[order(res_divided[,10]),]
  res_divided<-t(res_divided)
  
  library('pheatmap')
  pheatmap(res_divided,
           cluster_rows= F,cluster_cols= F,
           # color = colorRampPalette(c("white","#F5AC89","#BD2A34","#650520"))(1000),
           color = colorRampPalette(c("white","#9AC9E0","#2272B4","#1A325F"))(1000),
           border_color=NA #单元格边框设置
  )
  # 简化结果
  res_divided_three<-rbind(res_divided[1,],colSums(res_divided[2:9,]),res_divided[10,])
  pheatmap(res_divided_three,
           cluster_rows= F,cluster_cols= F,
           color = colorRampPalette(c("#1A325F","#2272B4","#9AC9E0","white","#F5AC89","#BD2A34","#650520"))(1000),
           border_color=NA #单元格边框设置
           )
  
}

# 获得circRNA的保守性得分
{
  filterExtremum<-function(x){
    # x <- c(23, 17, 14, 29, 37, 21, 16, 18, 26, 31)
    iqr <- IQR(x)
    upper_bound <- median(x) + 1.5 * iqr
    lower_bound <- median(x) - 1.5 * iqr
    x_filtered <- x[x >= lower_bound & x <= upper_bound]
    return(x_filtered)
  }
  
  df<-read.table('02Cos_GTF/04circ_diff/circ_all.pos.cos',header = F,stringsAsFactors = F,sep = '\t')
  colnames(df)<-c('type','chr','pos','value')
  df<-na.omit(df)
  df$type<-factor(df$type)
  
  library(dplyr)
  # 假设你的数据框为df，包含"group"和"value"两列
  result_mean <- df %>%
    group_by(type) %>%
    summarize(avg_value = mean(value))
  write.table(result_mean,'graph/cos_circ_average.txt',row.names = F,sep = '\t',quote = F)

  
  
  # 假设你的向量为vec
  res_list=list()
  for (i in unique(df$type)) {
    categories <- cut(df[df$type==i,'value'], breaks = seq(0, 1, 0.1), include.lowest = TRUE, labels = FALSE)
    res_list[[i]]=table(categories)
  }
  res<-do.call(rbind,res_list)
  # write.table(t(res),'graph/cos_circ_group.txt',sep = '\t',quote = F)
  res_divided <- apply(res,1, function(x) x/sum(x))
  # write.table(res_divided,'graph/cos_group_freq.txt',sep = '\t',quote = F)
  # 按照第10个模块的大小排序
  res_divided<-t(res_divided)
  res_divided<-res_divided[order(res_divided[,10]),]
  res_divided<-t(res_divided)
  
  library('pheatmap')
  pheatmap(res_divided,
           cluster_rows= F,cluster_cols= F,
           # color = colorRampPalette(c("white","#F5AC89","#BD2A34","#650520"))(1000),
           color = colorRampPalette(c("white","#9AC9E0","#2272B4","#1A325F"))(1000),
           border_color=NA #单元格边框设置
  )
  # 简化结果
  res_divided_three<-rbind(res_divided[1,],colSums(res_divided[2:9,]),res_divided[10,])
  pheatmap(res_divided_three,
           cluster_rows= F,cluster_cols= F,
           color = colorRampPalette(c("#1A325F","#2272B4","#9AC9E0","white","#F5AC89","#BD2A34","#650520"))(1000),
           border_color=NA #单元格边框设置
  )
}