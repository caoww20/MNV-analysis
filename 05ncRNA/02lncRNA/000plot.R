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
  mnv_diff<-fread('02RNAsnp/mnv_diff.res',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  snv_diff<-fread('02RNAsnp/snv_diff.res',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  mnv_snv_diff<-fread('02RNAsnp/merge_res.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
}
## 有多少MNV影响lncRNA的二级结构
{
  tmp<-mnv_diff[abs(mnv_diff$V4)==1,]
  tmp<-str_split_fixed(tmp$V1,'[|]',n=Inf)
  length(unique(tmp[,2]))
  # 总共有396MNV落在lncRNA上，标记为all
  length(unique(tmp[,2]))/426019
}
## 绘制数量的饼图
{
  # mnv
  a<-nrow(mnv_diff)
  # diff_ref
  b<-nrow(mnv_diff[abs(mnv_diff$V4)==1,])
  # diff_snv & diff_ref
  c<-nrow(mnv_snv_diff[mnv_snv_diff$V8==1,])
  
  # 745811 237008 75041
  
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
    filename = './graph/lncRNA_ratio1.pdf',
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
    filename = './graph/lncRNA_ratio2.pdf',
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
    filename = './graph/lncRNA_ratio.pdf',
    plot = p,width = 6,height = 6,
    units = 'in'
  )
}
## 绘制p的箱线图
{
  # install.packages("gghalves")
  library(gghalves)
  
  df <- data.frame(type = c(rep('mnv',nrow(mnv_diff)),rep('snv',nrow(snv_diff))), value = c(mnv_diff$V3,snv_diff$V3))
  df$type<-factor(df$type,levels = c('mnv','snv'))
  
  p<-ggplot(df, aes(x = type, y = value, fill = type))+
    geom_half_violin(aes(fill = type),
                     position = position_nudge(x = .15, y = 0),
                     adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
    geom_boxplot(aes(x = type,y = value, fill = type),
                 outlier.shape = NA,
                 width = .05,
                 color = "black")+
    labs(x=NULL,y='pvalue')+
    scale_fill_manual(values = c('#A53038','#3A5B9F'))+
    scale_y_continuous(limits = c(0,1))+
    theme_bw()+mythem
  
  # 保存
  ggsave(
    filename = './graph/lncRNA_p_change.pdf',
    plot = p,width = 6,height = 5,
    units = 'in'
  )
}
