# 载入包
library(data.table)
library(ggplot2)
library(ggforce)
# devtools::install_github('smin95/smplot2', force = TRUE)
library(smplot2)
## 主题
{
  mythem=theme(panel.grid=element_blank(),
               plot.margin=unit(rep(2,4),'lines'),
               legend.position="none",
               panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
               # linewidth 替换element_rect，element_rect已经弃用
               # title = element_text(vjust=1,size=16,face = "bold", color='black'),
               axis.text.x=element_text(vjust=1,size=16,face = "bold", colour='black'),
               axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
               axis.text.y=element_text(vjust=1,size=16,face = "bold", color = 'black'),
               axis.title.y=element_text(vjust=2, size=20,face = "bold", color='black')
  )
}
## 读入数据
{
  mnv_diff<-fread('02RNAFold/mnv/mnv_diff.res',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  snv_diff<-fread('02RNAFold/snv/snv_diff.res',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  mnv_snv_diff<-fread('02RNAFold/merge_res.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
}
## 有多少MNV影响pre-miRNA的二级结构
{
  tmp<-mnv_diff[abs(mnv_diff$V3)>1,]
  tmp<-str_split_fixed(tmp$V1,'[|]',n=Inf)
  length(unique(tmp[,2]))
  # 总共有396MNV落在pre-mIRNA上，标记为all
  length(unique(tmp[,2]))/396
}
## 绘制数量的饼图
{
  # mnv
  a<-nrow(mnv_diff)
  # diff_ref
  b<-nrow(mnv_diff[abs(mnv_diff$V3)>=1,])
  # diff_snv & diff_ref
  c<-nrow(mnv_snv_diff[abs(mnv_snv_diff$V3)>=1 & mnv_snv_diff$V5==1,])
  # 416 300 113
  
  # 大圆
    ratio <- c(1-b/a, b/a)
    type <- c('no_change', 'diff_ref')
    A <- data.frame(ratio, type)
    p<-ggplot()+
      xlab("")+ylab('')+
      geom_arc_bar(data=A,
                   stat = "pie",
                   aes(x0=0,y0=0,r0=1,r=2,
                       amount=ratio,fill=type)
      )+theme_bw()+mythem+scale_fill_manual(values = c( '#4396B1','#2A5888'))
    # 保存
    ggsave(
      filename = './graph/premiRNA_ratio1.pdf',
      plot = p,width = 6,height = 6,
      units = 'in'
    )
  # 小圆
    ratio <- c(1-c/b, c/b)
    type <- c('diff_ref','diff_snv')
    A <- data.frame(ratio, type)
    p<-ggplot()+
      xlab("")+ylab('')+
      geom_arc_bar(data=A,
                   stat = "pie",
                   aes(x0=0,y0=0,r0=1,r=2,
                       amount=ratio,fill=type)
      )+theme_bw()+mythem+scale_fill_manual(values = c( '#89CEED','#4396B1'))
    # 保存
    ggsave(
      filename = './graph/premiRNA_ratio2.pdf',
      plot = p,width = 6,height = 6,
      units = 'in'
    )
  # 合并圆
    ratio <- c(1-b/a, b/a-c/a, c/a)
    type <- c('no_change', 'diff_ref','diff_snv')
    A <- data.frame(ratio, type)
    p<-ggplot()+
      xlab("")+ylab('')+
      geom_arc_bar(data=A,
                   stat = "pie",
                   aes(x0=0,y0=0,r0=1,r=2,
                       amount=ratio,fill=type)
      )+theme_bw()+mythem+scale_fill_manual(values = c( '#89CEED','#4396B1','#2A5888'))
    # 保存
    ggsave(
      filename = './graph/premiRNA_ratio.pdf',
      plot = p,width = 6,height = 6,
      units = 'in'
    )

}
## 绘制G的箱线图
{
  df <- data.frame(type = c(rep('mnv',nrow(mnv_diff)),rep('snv',nrow(snv_diff))), value = c(mnv_diff$V3,snv_diff$V3))
  df$value<-abs(df$value)
  df$type<-factor(df$type,levels = c('mnv','snv'))
  
  p<-ggplot(df, aes(x = type, y = value, fill = type)) + 
    scale_y_continuous(breaks=seq(0, 15, by = 5),limits = c(0, 15))+
    labs(y="Absolute value of MFE change",x="")+
    scale_fill_manual(values = c('#A53038','#3A5B9F'))+
    sm_raincloud() +
    theme_bw()+mythem+
    scale_x_discrete(labels = c('MNV','SNV'))
  # 保存
  ggsave(
    filename = './graph/premiRNA_MFE_change.pdf',
    plot = p,width = 6,height = 5,
    units = 'in'
  )
}

