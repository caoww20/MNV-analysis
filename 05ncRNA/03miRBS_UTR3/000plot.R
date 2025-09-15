# 载入包
library(data.table)
library(ggplot2)
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

# UTR3
# 总共有多少MNV造成多少对差异的MNV-miRBS，gain多少，lose多少
{
  tmp<-fread('./UTR3/diff_MNV_vs_ref.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  tmp<-tmp[tmp$V3!='-' & tmp$V16!='unclear' & tmp$V15!='-',]
  mnv<-str_split_fixed(tmp$V1,'[|]',n=Inf)
  length(unique(mnv[,2]))
  tmp2<-cbind(mnv[,2],tmp[c(2,16)])
  tmp2<-tmp2[!duplicated(tmp2),]
  a<-tmp2
  remove(tmp,mnv,tmp2)
}
# 绘制条形图
{
  getNum<-function(file){
    data<-read.table(file,header = F,stringsAsFactors = F,sep = '\t')
    colnames(data)<-c('mnvid','gene')
    data<-data[!duplicated(data),]
    data<-as.data.frame(table(table(data$gene)),stringsAsFactors=F)
    colnames(data)<-c('id','number')
    data$id<-as.numeric(data$id)
    data[5,2]<-sum(data[data$id>=5,2])
    data<-data[data$id<=5,]
  }
  
  mnv<-getNum('./res/UTR3_MNVNum.txt')
  p<-ggplot(mnv,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#1A325F')+
    ylim(0, 6000)+
    labs(title='',x='Number of MNVs per Gene\'s UTR3',y='Number of Genes')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/UTR3_MNVNum.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  mnv_vs_ref<-getNum('./res/UTR3_MNVNum_vs_ref.txt')
  # 顺便查看大于等于5个MNV的数量
  p<-ggplot(mnv_vs_ref,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#2A5888')+
    ylim(0, 6000)+
    labs(title='',x='Number of miRBS changed MNVs per Gene\'s UTR3',y='Number of Genes')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/UTR3_MNVNum_vs_ref.pdf',
    plot = p,width = 6,height = 6,
    units = 'in'
  )
}
# 绘制大小饼图，有多少MNV-trans-miRBS区别于SNV-trans-miRBAS
# 如果是MNV-miRBS，则需要把trans去掉
{
  library(ggforce)
  # 绘制小饼图
  data<-fread('./UTR3/diff_MNV_vs_SNV.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  data<-table(data$V5)
  a<-data[1]+data[2]
  b<-data[1]
  ratio <- c(1-b/a, b/a)
  type <- c('no_change', 'change')
  A <- data.frame(ratio, type)
  
  p<-ggplot()+
    xlab("")+ylab('')+
    geom_arc_bar(data=A,
                 stat = "pie",
                 aes(x0=0,y0=0,r0=1,r=2, #r0=1表示空心环 
                     amount=ratio,fill=type, 
                     explode=c(-0.05)) # explode表示有空格
    )+theme_bw()+mythem+scale_fill_manual(values = c( '#4396B1','#2A5888'))
  ggsave(
    filename = './graph/UTR3_ratio.pdf',
    plot = p,width = 6,height = 6,
    units = 'in'
  )
}
## 冲击图
{
  library(alluvial)
  data<-fread('./res/UTR3_alluvial.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
  data<-data[data$Freq>0,]
  colnames(data)<-c('snv','mnv','freq')
  # 绘图
  rownames(data)<-NULL
  data<-data[order(data$snv),]
  write.table(data,'./graph/UTR3_alluvial_2.txt',quote = F,sep = '\t',row.names = F)
  
  set.seed(8) # for nice colors
  # cols <- hsv(h = sample(1:8/10), s = sample(3:8)/8, v = sample(3:8)/8)
  ord <- list(NULL, with(data, order(snv, mnv)))
  alluvial(
    data[,1:2], alpha=1,
    freq = data$freq,
    blocks = T,
    gap.width=0.1, # 块的间隔
    xw = 0.3, # 线的曲度
    cw = 0.01, # 字块大小
    # col=c('#BDBDBD','#4EAE4A','#984EA3'),
    # axis_labels='',
    axes=F,
    ann=F,
    # col = cols[match(data$snv, unique(data$mnv))],
    ordering=ord
  )
  
  #600*600
}

# lncRNA
# 总共有多少MNV造成多少对差异的MNV-miRBS，gain多少，lose多少
{
  tmp<-fread('./lncRNA/diff_MNV_vs_ref.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  tmp<-tmp[tmp$V3!='-' & tmp$V16!='unclear' & tmp$V15!='-',]
  mnv<-str_split_fixed(tmp$V1,'[|]',n=Inf)
  length(unique(mnv[,2]))
  tmp2<-cbind(mnv[,2],tmp[c(2,16)])
  tmp2<-tmp2[!duplicated(tmp2),]
  b<-tmp2
  remove(tmp,mnv,tmp2)
}
# 绘制条形图
{
  mnv<-getNum('./res/lncRNA_MNVNum.txt')
  p<-ggplot(mnv,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#4396B1')+
    ylim(0, 25000)+
    labs(title='',x='Number of MNVs per lncRNA',y='Number of lncRNAs')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/lncRNA_MNVNum.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  mnv_vs_ref<-getNum('./res/lncRNA_MNVNum_vs_ref.txt')
  p<-ggplot(mnv_vs_ref,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#89CEED')+
    ylim(0, 25000)+
    labs(title='',x='Number of miRBS changed MNVs per lncRNA',y='Number of lncRNAs')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/lncRNA_MNVNum_vs_ref.pdf',
    plot = p,width = 6,height = 6,
    units = 'in'
  )
}
# 绘制饼图，有多少MNV-trans-miRBS区别于SNV-trans-miRBAS
{
  data<-fread('./lncRNA/diff_MNV_vs_SNV.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  data<-table(data$V5)
  a<-data[1]+data[2]
  b<-data[1]
  library(ggforce)
   <- c(1-b/a, b/a)
  type <- c('no_change', 'change')
  A <- data.frame(ratio, type)
  
  p<-ggplot()+
    xlab("")+ylab('')+
    geom_arc_bar(data=A,
                 stat = "pie",
                 aes(x0=0,y0=0,r0=1,r=2, #r0=1表示空心环 
                     amount=ratio,fill=type, 
                     explode=c(-0.05)) # explode表示有空格
    )+theme_bw()+mythem+scale_fill_manual(values = c( '#89CEED','#4396B1'))
  ggsave(
    filename = './graph/lncRNA_ratio.pdf',
    plot = p,width = 6,height = 6,
    units = 'in'
  )
}
## 冲击图
{
  library(alluvial)
  data<-fread('./res/lncRNA_alluvial.txt',header = F,stringsAsFactors = F,sep = '\t',data.table = F)
  data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
  data<-data[data$Freq>0,]
  colnames(data)<-c('snv','mnv','freq')
  # 绘图
  rownames(data)<-NULL
  data<-data[order(data$snv),]
  write.table(data,'./graph/lncRNA_alluvial_2.txt',quote = F,sep = '\t',row.names = F)
  
  set.seed(8) # for nice colors
  # cols <- hsv(h = sample(1:8/10), s = sample(3:8)/8, v = sample(3:8)/8)
  ord <- list(NULL, with(data, order(snv, mnv)))
  alluvial(
    data[,1:2], alpha=1,
    freq = data$freq,
    blocks = T,
    gap.width=0.1, # 块的间隔
    xw = 0.3, # 线的曲度
    cw = 0.01, # 字块大小
    # col=c('#BDBDBD','#4EAE4A','#984EA3'),
    # axis_labels='',
    axes=F,
    ann=F,
    # col = cols[match(data$snv, unique(data$mnv))],
    ordering=ord
  )
  
  #600*600
}

# 所有的miRBS的结果
c<-rbind(a,b)
c<-c[!duplicated(c),]

