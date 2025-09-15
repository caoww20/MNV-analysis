## 主题
{
  library(ggplot2)
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
# remotes::install_github("mbojan/alluvial", build_vignettes=TRUE)
library(alluvial)
library(data.table)


## 绘制alluvial
{
  # 读入数据
  data<-read.table('02canonicalTrans/alluvial.txt',header = F,sep = '\t',stringsAsFactors = F)
  data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
  data<-data[data$Freq>0,]
  colnames(data)<-c('snv','mnv','freq')
  
  # 修饰数据（但龚老师不同意）
  {
    # 考虑到有两个值特别大，所以进行适当的缩放
    data[1,3]=3000
    data[3,3]=2500
    data[4,3]=2500
    # data<-data[-10,] # 对于synonymous_variant->stop_lost,新的版本已经归到stop_lost->stop_lost，所以这里就不需要了
    data[data$freq<=48,3]<-100
  }

  # 绘图
  rownames(data)<-NULL
  data<-data[order(data$snv),]
  data[1,2]<-'c_missense_variant_1'
  data[2,2]<-'c_missense_variant_2'
  
  set.seed(8) # for nice colors
  cols <- hsv(h = sample(1:8/10), s = sample(3:8)/8, v = sample(3:8)/8)
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
    col = cols[match(data$snv, unique(data$mnv))],
    ordering=ord
  )
  #800*400
}

# AA
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
  
  AA<-getNum('03AA/AA_perGene.txt')
  p<-ggplot(AA,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#1A325F')+
    ylim(0, 5000)+
    labs(title='',x='Number of amino acid altered MNVs per Gene',y='Number of Genes')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/AA.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  stop_change<-getNum('03AA/stop_change_perGene.txt')
  p<-ggplot(stop_change,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#2A5888')+
    ylim(0, 1600)+
    labs(title='',x='Number of stop changed MNVs per Gene',y='Number of Genes')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/stop_change.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  star_change<-getNum('03AA/start_change_perGene.txt')
  p<-ggplot(star_change,aes(x = id,y = number))+
    geom_col(position = 'dodge',width = 0.8,fill='#4396B1')+
    ylim(0, 1200)+
    labs(title='',x='Number of start change MNVs per Gene',y='Number of Genes')+
    theme_bw()+mythem
  ggsave(
    filename = 'graph/start_change.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}

## PTV的数量
{
  data<-read.table('02canonicalTrans/merge_total.txt',header = T,sep = '\t',stringsAsFactors = F)
  data<-data[c(3,15,19,23,27:37)]
  # complete gain stop。两个点都没有获得，但是MNV获得终止
  sum(data$complete.gain.stop) # 522
  # rescue stop。只有是拯救的都是
  sum(data$rescue.stop) # 1986
}









