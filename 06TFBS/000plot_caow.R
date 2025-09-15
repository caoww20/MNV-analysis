library(tidyverse)
# 主题
{
  mythem=theme(panel.grid=element_blank(),
               plot.margin=unit(rep(2,4),'lines'),
               legend.position="none",
               panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
               # linewidth 替换element_rect，element_rect已经弃用
               # title = element_text(vjust=1,size=16,face = "bold", color='black'),
               axis.text.x=element_text(vjust=1,size=8,face = "bold", colour='black'),
               axis.title.x=element_text(vjust=0, size=8,face = "bold", color='black'),
               axis.text.y=element_text(vjust=0.5,size=8,face = "bold", color = 'black'),
               axis.title.y=element_text(vjust=2, size=8,face = "bold", color='black')
  )
}
# 2个库的TF重叠情况 venn
{
  library(eulerr)
  fit1 <- euler(c("A" = 707, "B" = 401, 
                  "A&B" = 347))
  
  fill1=list(fill=c("#e6194B","#5B9BD5"),ncol=2,
             alpha=1)
  
  plot(fit1,fills = fill1,
       #labels = list(col="white",font=3,cex=2),
       quantities = list(col="black",cex=2),
       legend = list(labels=c("JASPAR","HOCOMOCO"),cex=2,col="white",fontsize=10,side='top',pch=21),
       edges = list(lwd=3,lty=2)
  )
}
# 计算两个库的motif的相似度，相似度90%以上的有多少，p值为0.01
{
  motif_sim = read_tsv("motif_sim",comment = '#') # motif相似性表
}
# 转录本启动子区域平均有多少MNV（unique MNVs/去冗余后的up长度），分别计算所有转录本、protein_trans(classic)的和lncRNA_trans（选用MNV的AC大于1）
{
  
}
# 每个转录本启动子区域携带MNV的数量分布，如携带1个MNV的转录本有多少，2个的有多少...（选用MNV的AC大于1，这里指的是转录本）
# 使用的是经典转录本
{
  mnv_promoter=read_tsv('mnv_fiilter_ac1_p_promoter',col_names = F)
  # 统计十个以内的MNV
  mnvl10=mnv_promoter%>% group_by(X8)%>% summarise(n=n()) %>% group_by(n) %>% summarise(n1=n()) %>% filter(n<10)
  mnvm10=mnv_promoter%>% group_by(X8)%>% summarise(n=n()) %>% group_by(n) %>% summarise(n1=n()) %>% filter(n>9)
  sum(mnvm10$n1) # 36002
  rbind(mnvl10,c('>10',sum(mnvm10$n1))) %>%
    mutate(n1=as.integer(n1)) %>% mutate(n=ordered(n,levels=n)) %>% 
    ggplot(aes(n,n1))+
    geom_col(position = 'dodge',width = 0.8, fill='#77398E')+
    ylim(0, 80000)+
    labs(title='',x='Number of MNVs',y='Number of transcripts')+
    theme_bw()+mythem
}
# 合并两个库的结果，直接合并，因为他们的motif不同，所以不进行去重处理
{
  # 载入JASPAR HOCOMOCO的结果
  x2 = read_tsv("jaspar_hg38",col_types = cols(Positions = col_character()))
  x3 = read_tsv("hocomoco_hg38",col_types = cols(Positions = col_character()))
  mnvid=read_tsv('mnvlist_hg38_final_version',col_names = F,col_types = cols(X2 = col_character()))
  colnames(mnvid)=c('Chr', 'Positions','MNVID','Refs',  'Alts',  'rsIDs')
  x2 %>% inner_join(mnvid,by=c('Chr', 'Positions','Refs',  'Alts',  'rsIDs'))
  all_result = x2 %>% bind_rows(x3)
  all_result=all_result %>% inner_join(mnvid,by=c('Chr', 'Positions','Refs',  'Alts',  'rsIDs')) %>% mutate(MNVID.x=MNVID.y) %>% select(-MNVID.y)
  colnames(all_result)[5]='MNVID'
  write.table(all_result,file = 'all_result_h_j',quote = F,col.names = T,row.names = F,sep = "\t")
}
# FIMO鉴定的结果统计，两个方法两种改变 （gain、loss的条形图）
{
  library(reshape2)
  gainLoss<-list()
  gainLoss[['hocomoco']]<-table(x3$Effect)
  gainLoss[['jaspar']]<-table(x2$Effect)
  gainLoss<-as.data.frame(do.call(rbind,gainLoss))
  gainLoss$type<-rownames(gainLoss)
  gainLoss<-melt(gainLoss,id='type')
  colnames(gainLoss)[2]='Effect'
  
  ggplot(gainLoss, aes(x = type, y = value, fill = Effect)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Method") +
    ylab("Number of MNV-TFBS pairs") +
    theme_bw() + mythem+
    theme(legend.position='top',
          legend.key.size = unit(30,'pt'))+  theme(legend.text = element_text(size=16,face = "bold", colour='black'))+  # 调整图例的字体大小
    scale_fill_manual(values = c("#1A325F","#ac2e31"))
}
# 总共有多少个MNV（xx%）造成了多少个MNV-TFBS的改变
{
  
}
# 总共有多少个MNV-TFBS（xx%）造成与SNV-TFBS不同的改变
{
  
}
# 绘制冲击图 表示MNV和SNV对TFBS效应的差异
{
  # colnames(data)<-c('snv','mnv')
  # table(data) # 为3*3组合数
  library(alluvial)
  library(data.table)
  data<-fread('snv_mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)
  data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
  data<-data[data$Freq>0,]
  colnames(data)<-c('snv','mnv','freq')
  data<-data[order(data$snv),]
  
  ord <- list(NULL, with(data, order(snv, mnv)))
  alluvial(
    data[,1:2], alpha=0.95,
    freq = data$freq,
    blocks = T,
    gap.width=0.1, # 块的间隔
    xw = 0.3, # 线的曲度
    cw = 0.01, # 字块大小
    axes=F,
    ann=F,
    col = c("#469566","#192d5a","#754f84","#192d5a","#754f84", "#469566", "#192d5a"),
    ordering=ord
  )
  
}
# 正负2.5kb分块，100bp为一单位(使用的jaspar？用的是所有转录本的数据？)
{
  # MNV-TFBS的数量曲线
  # 差异的MNV-TFBS的数量曲线
  # CPG的曲线
  
  # MNV-TFBS的数量曲线
  {
    jaspar_can = read_tsv('jaspar_5k_promoter_hg38',col_types = cols(Positions = col_character()))
    mnv_dis_p_c = read_tsv('jaspar_mnv_distance_promoter_5k',col_names = F)
    colnames(mnv_dis_p_c)=c('MNVID','gene','dis')
    
    dataplot1=jaspar_can %>% select(TF,MNVID) %>% inner_join(mnv_dis_p_c,by='MNVID') %>% filter(dis>0) %>% group_by(TF,MNVID) %>% 
      summarise(n=n(),dis=mean(dis)) %>% filter(n==1) %>%  mutate(a=cut(dis,breaks = seq(0,5000,by=100))) %>% group_by(a) %>% 
      summarise(num=n()) %>% filter(!is.na(a))
    
    # 统计MNV-TFBSpairs的数量
    point2 = spline(dataplot1$a,dataplot1$num,5000)
    abab=data.frame(point2$x,point2$y)
    x=abab$point2.x
    new_x <- scale(x, center = min(x), scale = max(x) - min(x)) * (2500 - (-2500)) + (-2500)
    
    abab$point2.x=new_x
    p1 = ggplot(abab,aes(x =  point2.x)) +
      geom_point(aes(y = point2.y),col='#BE1D2C')+
      scale_x_continuous(breaks=seq(-2500,2500,500))+
      ylim(c(0,30000))+theme_bw()+mythem+theme(plot.margin=unit(rep(1,4),'lines'))
  }
  # 统计启动子区域MNV影响的TFBS数量
  {
    mnvdis=mnv_dis_p_c%>% group_by(MNVID) %>%  summarise(n=n(),dis=mean(dis)) %>% filter(n==1)%>%select(MNVID,dis)
    # tfbs_permnv  = jaspar_can %>% select(TF,MNVID)%>%group_by(MNVID)%>% count() %>%
    #   inner_join(mnvdis,by='MNVID')%>%  mutate(a=cut(dis,breaks = seq(0,5000,by=100)))
    tfbs_permnv  = jaspar_can %>% select(TF,MNVID)%>%group_by(MNVID)%>% summarise(n=n()) %>%
      inner_join(mnvdis,by='MNVID')%>%  mutate(a=cut(dis,breaks = seq(0,5000,by=500)))
    
    p2 = ggplot(tfbs_permnv,aes(color = a))+ 
      geom_boxplot(aes(x = a, y =n),col='#177e7c',outlier.shape = NA,position = position_nudge(x = 0.5)
      )+ coord_cartesian(ylim = c(0, 35))+scale_x_discrete(labels=seq(-2500,2600,500))+
      theme_bw()+mythem+theme(plot.margin=unit(rep(1,4),'lines')) # 设置X轴标签的方向
    
  }

  # 合并结果
  {
    library(gridExtra)
    g1 <- ggplotGrob(p1)
    g2 <- ggplotGrob(p2)
    g3 <- ggplotGrob(p3)
    
    grid.arrange(rbind(g3, g1, g2, size = "last"))
  }
}







