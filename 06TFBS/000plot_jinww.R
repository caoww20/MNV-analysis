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
  get_gradient_colors<-function(colors,n){
    # colors <- c("#FF0000", "#0000FF")
    
    # 生成渐变色函数
    gradient_colors <- colorRampPalette(colors)
    
    # 生成10个渐变色
    # n <- 10
    result <- gradient_colors(n)
    return(result)
  }
}
# 查看两个数据库motif的相似性
{
  df<-read.table('old/motif_sim',header = T,stringsAsFactors = F,sep = '\t')
  motif<-read.table('../01motif/jaspar_motif_id',header = F,sep = ' ',stringsAsFactors = F)
  df$diff<-df$Query_consensus==df$Target_consensus
  res<-merge(df,motif[2:3],by.x='Query_ID',by.y='V2',sort = F,all.x = T)
  
}
# 查看lncRNA，mRNA，all的MNV-pair数，注意位置的区别可能导致同时有获得和丢失
{
  lnc<-fread('./04res_lnc/jaspar_pair.txt',header = F,stringsAsFactors = F,sep='\t',data.table = F)
  mRNA<-fread('./04res_mRNA/jaspar_pair.txt',header = F,stringsAsFactors = F,sep='\t',data.table = F)
  all<-rbind(lnc,mRNA)
  all<-all[!duplicated(all),]
  c(nrow(lnc),nrow(mRNA),nrow(all)) #  2440476 790532 3123623
  remove(lnc,mRNA,all)
}
# 统计发现，编码基因启动子区域中平均含有xx个MNV，其中xx%含有至少5个及其以上的MNV改变。
{
  df<-read.table('./04res_mRNA/geneMNVNum.txt',header = F,sep = '\t',stringsAsFactors = F)
  mean(df$V2) # 12
  median(df$V2) # 8
  nrow(df[df$V2>=5,])/nrow(df) # 72.25%
}
# 平均/中位数每个基因含有xx个MNV影响TFBS结合，xx个cahnged MNV-TFBS pair
{
  a<-read.table('./04res_mRNA/jaspar_gene.num',header = F,sep = '\t',stringsAsFactors = F)
  # 含有xx个MNV影响TFBS结
  mean(a$V3) # 9
  # 含有差异的MNV-TFBS pair的数量
  mean(a$V4) # 39
  # 绘制条形图
  {

    data<-as.data.frame(table(a$V3),stringsAsFactors=F)
    colnames(data)<-c('id','number')
    data$id<-as.numeric(data$id)
    data<-data[2:nrow(data),]
    data[10,2]<-sum(data[data$id>=10,2])
    data<-data[data$id<=10,]
    data$id<-factor(data$id)
  
    p<-ggplot(data,aes(x = id,y = number))+
      geom_col(position = 'dodge',width = 0.8,fill='#1A325F')+
      ylim(0, 6000)+
      labs(title='',x='Number of MNVs per gene\'s promotor',y='Number of genes')+
      theme_bw()+mythem
    ggsave(
      filename = 'graph/tfbs-altered_MNV.pdf',
      plot = p,width = 8,height = 6,
      units = 'in'
    )
      }
  }
# 2个方法得到的gain loss结果
{
  library(reshape2)
  colnames(df)<-c('gain','loss')
  df<-reshape2::melt(df,by=0)
  colnames(df)<-c('type','Effect','value')
  
  p<-ggplot(df, aes(x = type, y = value, fill = Effect)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("The Library of TF Motifs") +
    ylab("Number of MNV-TFBS pairs") +
    scale_y_continuous(breaks = seq(0, 450000, 150000), limits = c(0, 450000)) +
    theme_bw() + mythem+
    scale_fill_manual(values = c("#2A5888","#4396B1"))
  ggsave(
    filename = './graph/mRNA_tfbs.pdf',
    plot = p,width = 6,height = 5,
    units = 'in'
  )
}
# 我们也进一步评估了区别于SNV的影响的MNV，发现xx%的MNV-TFBS区别于SNV-TFBS pair，也通过冲击图进一步展示了这种变化。
{
  library(alluvial)
  library(data.table)
  {
    data<-fread('./04res_mRNA/jaspar_mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)
    table(data$V2)[2]/nrow(data) 
    a<-fread('./04res_mRNA/jaspar_mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)
  }
  # 绘制冲击图
  {
    data<-fread('./04res_mRNA/jaspar_snv_mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)
    data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
    data<-data[data$Freq>0,]
    colnames(data)<-c('snv','mnv','freq')
    data<-data[order(data$snv),]
    write.table(data,'./graph/mRNA_jaspar_alluvial.txt',sep = '\t',quote = F,row.names = F)
        
    data<-fread('./04res_mRNA/jaspar_snv_mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)
    data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
    data<-data[data$Freq>0,]
    colnames(data)<-c('snv','mnv','freq')
    data<-data[order(data$snv),]
    write.table(data,'./graph/mRNA_merge_alluvial.txt',sep = '\t',quote = F,row.names = F)
    
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
      # col = c("#469566","#192d5a","#754f84","#192d5a","#754f84", "#469566", "#192d5a"),
      ordering=ord
    )
  }
}

# 启动子区域CPG island的分布
{
  cpgi<-read.table('./03caow/cpgi_can_promoter.txt',sep='\t',stringsAsFactors = F)
  colnames(cpgi)=c('f')
  dist=seq(-2500,2499)
  cpgi=data.frame(cpgi$f,dist)
  p=ggplot(cpgi,aes(x =  dist)) +
    geom_point(aes(y = cpgi.f),col='#425190')+
    scale_x_continuous(breaks=seq(-2500,2500,500))+
    ylim(c(0,1))+theme_bw()+mythem
  ggsave(
    filename = './graph/mRNA_cpgi_per_region.pdf',
    plot = p,width = 8,height = 4,
    units = 'in'
  )
}

# 统计启动子区域MNV影响的TFBS数量(经典转录本)
{
  jaspar_can = read_tsv('./03caow/jaspar_5k_promoter_hg38',col_types = cols(Positions = col_character()))
  
  mnv_dis_p_c_j = read_tsv('./03caow/jaspar_mnv_distance_promoter_5k',col_names = F)
  colnames(mnv_dis_p_c_j)=c('MNVID','gene','dis')
  mnv_dis_p_c_h = read_tsv('./03caow/hocomoco/hocomoco_mnv_distance_promoter_5k',col_names = F)
  colnames(mnv_dis_p_c_h)=c('MNVID','gene','dis')
  # 统计启动子区域MNV影响的TFBS数量
  mnvdis_j=mnv_dis_p_c_j%>% group_by(MNVID) %>%  summarise(n=n(),dis=mean(dis)) %>% filter(n==1)%>%select(MNVID,dis)
  mnvdis_h=mnv_dis_p_c_h%>% group_by(MNVID) %>%  summarise(n=n(),dis=mean(dis)) %>% filter(n==1)%>%select(MNVID,dis)
  
  tfbs_permnv_j  = jaspar_can %>% select(TF,MNVID)%>%group_by(MNVID)%>%count() %>%
    inner_join(mnvdis_j,by='MNVID')%>%  mutate(a=cut(dis,breaks = seq(0,5000,by=500)))

  p2<-ggplot(tfbs_permnv_j,aes(color = a))+ 
    geom_boxplot(aes(x = a, y =n),col='#177e7c',outlier.shape = NA,position = position_nudge(x = 0.5)
    )+ coord_cartesian(ylim = c(0, 35))+
    scale_x_discrete(labels=seq(-2500,2600,500))+
    theme_bw()+mythem 
  
  library(gghalves)
  p2<-ggplot(tfbs_permnv_j, aes(x = a, y = n, fill = a))+
    geom_half_violin(aes(fill = a),
                     position = position_nudge(x = .15, y = 0),
                     adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
    geom_boxplot(aes(x = a,y = n, fill = a),
                 outlier.shape = NA,
                 width = .05,
                 color = "black")+
    labs(x=NULL,y='Affected TFBSs per MNV')+
    scale_fill_manual(values = get_gradient_colors(c('#2A5888','#C4F5FC'),10) )+
    scale_x_discrete(labels=seq(-2500,2600,500))+
    scale_y_continuous(limits = c(0,35))+
    theme_bw()+mythem
  ggsave(
    filename = './graph/mRNA_tfbs_per_MNV.pdf',
    plot = p2,width = 8,height = 4,
    units = 'in'
  )
}

