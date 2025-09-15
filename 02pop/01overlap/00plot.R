library(ggplot2)
library(ggsci)
library(scales)
library(data.table)
#提取AAAS配色
{
  mythem<-theme(panel.grid=element_blank(),
                plot.margin=unit(rep(2,4),'lines'),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                title = element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.text.x=element_text(vjust=1,size=16,angle = 0, face = "bold", color='black'),
                axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black'),
                legend.position="none"
                # legend.title = element_blank()
  )
  mythem_legend<-theme(panel.grid=element_blank(),
                       plot.margin=unit(rep(2,4),'lines'),
                       panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                       title = element_text(vjust=1,size=16,face = "bold", color='black'),
                       axis.text.x=element_text(vjust=1,size=16,angle = 0, face = "bold", color='black'),
                       axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                       axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                       axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black')
                       # legend.position="none"
                       # legend.title = element_blank()
  )
}

# 计算不同人群的MNV分组情况
{
  a<-read.csv('venn/venn_1000G.txt',header = F,stringsAsFactors = F,sep = '\t')
  x <- list(AFR=a[a$V1=='AFR',2],
            AMR=a[a$V1=='AMR',2],
            EAS=a[a$V1=='EAS',2],
            EUR=a[a$V1=='EUR',2],
            SAS=a[a$V1=='SAS',2])
  venn::venn(x,zcolor=2:5,box=T)
  # 手动保存为图1000G_overlap.pdf 8*8
  table(a$V1)
  # AFR     AMR     EAS     EUR     SAS 
  # 1840367 1203314 888644 979825 1048132 
  # uniq
  # 808578 178624 263455 167428 287779
  # 43.93% 14.84% 29.65% 17.09% 27.46%
}

#########因为上面得到的uniq很多，所以我们进一步看了SNV的情况和MNV uniqu的分布###########################################

# 计算不同人群的SNV分组情况
{
  a<-read.table('snv/snv_overlap.txt',header = T,sep = ' ',stringsAsFactors = F)
  a<-a[order(-a$X1000G),]
  a$flag<-gsub('.snv','',a$pop)
  b<-reshape2::melt(a[,c('unique','notUnique','flag')],id='flag')
  b$flag<-factor(b$flag, levels = c('AFR','AMR','EAS','EUR','SAS'))
  
  p<-ggplot(data=b,aes(x=flag,y=value,fill=variable))+
    geom_bar(stat = 'identity',show.legend = F)+
    labs(title='',x='',y='')+
    scale_fill_manual(values = c("unique" = "#B0282E", "notUnique" = "#2A5888"))+
    scale_y_continuous(breaks = seq(0, 40000000, 10000000), limits = c(0, 40000000), expand = c(0,0))+
    theme_bw()+mythem
  ggsave(
    filename = '../graph/1000G_snv_overlap.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  # pop percent
  # 1 AFR.snv     39%
  # 2 AMR.snv     13%
  # 3 SAS.snv     29%
  # 4 EUR.snv     15%
  # 5 EAS.snv     33%
}

# 计算不同人群的common_SNV涉及到的MNV分组情况
{
  a<-read.csv('venn/venn_1000G_common.txt',header = F,stringsAsFactors = F,sep = '\t')
  x <- list(AFR=a[a$V1=='AFR',2],
            AMR=a[a$V1=='AMR',2],
            EAS=a[a$V1=='EAS',2],
            EUR=a[a$V1=='EUR',2],
            SAS=a[a$V1=='SAS',2])
  venn::venn(x,zcolor=2:5,box=T)
  # 手动保存到1000G_overlap_common.pdf 8*8
  table(a$V1)
  # AFR    AMR    EAS    EUR    SAS 
  # 459661 426229 427428 430886 435297 
  # length(unique(a$V2)) 541397
  # uniq
  # 30746 10138 15515 13136 15152
  # 6.69% 2.38% 3.63% 3.05% 3.48%
}


