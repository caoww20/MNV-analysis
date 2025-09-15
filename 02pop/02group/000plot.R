library(ggplot2)
# BiocManager::install("qqman") #利用BiocManager安装
library("qqman") #加载R包
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
# pop的信息
{
  pop<-read.table('/data/jinww/04reference/publicDB/1000G/igsr-1000_genomes_30x_on_grch38.tsv',header = T,sep = '\t')
  k<-as.data.frame(table(pop[,c(4,6)]))
  k<-k[k$Freq>2,]
  # pop<-pop[!duplicated(pop),]
  # colnames(pop)<-c('V1','V2')
  pop<-pop[c(1,4,6)]
  n=pop[pop$Superpopulation.code!='AMR',c(1,3)]
  m=pop[pop$Superpopulation.code=='AMR',c(1,2)]
  colnames(n)<-c('V2','V1')
  colnames(m)<-c('V2','V1')
  pop<-rbind(n,m)
}
# SNV和MNV的PCA图片
{
  # MNV PCA
  df<-read.table('04pca/pca_MNVFromCommonSNVFilter.eigenvec',header = F,stringsAsFactors = F)
  df<-df[2:4]
  pop<-read.table('01origin/sample/pop2.sample',header = F,stringsAsFactors = F)
  res<-merge(df,pop,by.x='V2',by.y='V2',all.x=T)
  res<-na.omit(res)
  colnames(res)<-c('sample','pca1','pca2','pop')
  p<-ggplot(res, aes(x=pca1, y=pca2, color=pop)) +
    geom_point(size = 3, alpha = 0.5)+
    scale_color_manual(values = c("#7E277B", "#FCA41C","#0D8040","#3353A3","#EF251F"))+
    theme_bw()+mythem_legend
  ggsave(
    filename = '../graph/pca_MNV.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  # scale_color_manual(values = c("#7E277B", "#FCA41C","#0D8040","#3353A3","#EF251F","#BA54A2"))
  # SNV PCA
  df<-read.table('04pca/pca_1000G.eigenvec',header = F,stringsAsFactors = F)
  df<-df[2:4]
  pop<-read.table('01origin/sample/pop2.sample',header = F,stringsAsFactors = F)
  res<-merge(df,pop,by.x='V2',by.y='V2',all.x=T)
  res<-na.omit(res)
  colnames(res)<-c('sample','pca1','pca2','pop')
  p<-ggplot(res, aes(x=pca1, y=pca2, color=pop)) +
    geom_point(size = 3, alpha = 0.5)+
    scale_color_manual(values = c("#7E277B", "#FCA41C","#0D8040","#3353A3","#EF251F"))+
    theme_bw()+mythem_legend
  ggsave(
    filename = '../graph/pca_SNV.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
# 计算SNV和MNV的FST区间曼哈顿图
{
  # MNV FST
  df<-read.table('03fst/filter/all_100k_5.windowed.weir.fst',header = T,stringsAsFactors = F)
  df<-cbind(paste0('snp',rownames(df)),df[c(1,2,5)])
  colnames(df)<-c('SNP','CHR','POS','Fst')
  manhattan(df,chr="CHR",bp="POS",p="Fst",snp="SNP", col=c("#3353A3", "#6CC6D6"),chrlabs = as.character(c(1:22)),logp=F,suggestiveline = F, 
    genomewideline = F,ylab="Fst",ylim=c(0,1),font.lab=4,cex.lab=1.2,main="MNV FST",cex=0.8)
  # SNV FST
  df<-read.table('03fst/snv/all_100k.windowed.weir.fst',header = T,stringsAsFactors = F)
  df<-cbind(paste0('snp',rownames(df)),df[c(1,2,5)])
  colnames(df)<-c('SNP','CHR','POS','Fst')
  df$CHR<-as.numeric(gsub('chr','',df$CHR))
  manhattan(df,chr="CHR",bp="POS",p="Fst",snp="SNP", col=c("#A71F36", "#F58E7E"),chrlabs = as.character(c(1:22)),logp=F,suggestiveline = F, 
            genomewideline = F,ylab="Fst",ylim=c(0,1),font.lab=4,cex.lab=1.2,main="SNV FST",cex=0.8)
}
# 平均权重
{
  # MNV FST
  df<-read.table('03fst/filter/all_100k_5.windowed.weir.fst',header = T,stringsAsFactors = F)
  # 0.686273
  df<-read.table('03fst/filter/all_5.weir.fst',header = T,stringsAsFactors = F)
  # 0.755189
  
  # SNV FST
  df<-read.table('03fst/snv/all_100k.windowed.weir.fst',header = T,stringsAsFactors = F)
  # 0.532347
  df<-read.table('03fst/snv/all.weir.fst',header = T,stringsAsFactors = F)
  # 0.920606
}
# 计算两者的相关性
{
  mnv<-read.table('03fst/filter/all_100k_5.windowed.weir.fst',header = T,stringsAsFactors = F)
  mnv<-mnv[c(1,2,3,5)]
  mnv$id<-paste0(mnv$CHROM,'_',mnv$BIN_START)
  snv<-read.table('03fst/snv/all_100k.windowed.weir.fst',header = T,stringsAsFactors = F)
  snv$CHROM<-gsub('chr','',snv$CHROM)
  snv<-snv[c(1,2,3,5)]
  snv$id<-paste0(snv$CHROM,'_',snv$BIN_START)
  res<-merge(snv,mnv[4:5],by='id',all.x=T,sort=F)
  res<-na.omit(res)
   
  dfcor<-cor.test(res$WEIGHTED_FST.x,res$WEIGHTED_FST.y,method = 'spearm') # rho=0.6594451
  # print(c(dfcor$estimate,dfcor$p.value))
  # ggplot(data=res,aes(x=WEIGHTED_FST.x,y=WEIGHTED_FST.y))+
  #   geom_point(pch=16,size=3,alpha=0.5)+
  #   labs(title = '',x='',y='')+ 
  #   theme_bw() + mythem
}
# 绘制MNV和SNV人群的热图以及相减图
{
  library(pheatmap)
  mnvHeatmap = read.table("./03fst/significant/heatmap_mnv.txt",header = T,row.names = 1)
  # my_palette <- colorRampPalette(c("blue", "red"))(n = 100)
  # my_colors = colorRampPalette(c("#1A325F","#2272B4","#9AC9E0","white","#F5AC89","#BD2A34","#650520"))(100)
  # my_colors <- colorRampPalette(c("#9AC9E0","#2272B4","#1A325F","#040419"))(100)
  # my_colors <- colorRampPalette(colors = c("#F06143","#E33641","#AE1758","#471D49","#040419"))(100)
  my_colors <- colorRampPalette(c("#9AC9E0","#2272B4","#1A325F"))(100)
  my_colors <- colorRampPalette(colors = c("#F06143","#E33641","#AE1758","#471D49"))(100)
  pheatmap(mnvHeatmap,
           display_numbers = TRUE,
           number_color = "white",
           number_format = "%.3f",
           border_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           fontsize=16,
           color = my_colors,
           breaks = seq(0, 0.2, by = 0.002) # breaks的数量最好和my_colors一致
           )
  snvHeatmap = read.table("./03fst/significant/heatmap_snv.txt",header = T,row.names = 1)
  pheatmap(snvHeatmap,
           display_numbers = TRUE,
           number_color = "white",
           number_format = "%.3f",
           border_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           fontsize=16,
           color = my_colors,
           breaks = seq(0, 0.2, by = 0.002) # breaks的数量最好和my_colors一致
  )
  # diffHeatmap = read.table("./03fst/significant/heatmap_diff.txt",header = T,row.names = 1)
  # pheatmap(diffHeatmap)
  # pheatmap(snvHeatmap,
  #          display_numbers = TRUE,
  #          number_color = "white",
  #          number_format = "%.3f",
  #          border_color = "white",
  #          fontsize=16,
  #          color = my_colors,
  #          breaks = seq(0, 0.2, by = 0.002) # breaks的数量最好和my_colors一致
  # )
}
# 计算不同人群的FST≥0.15的MNV的数量（以及任一都大于0.15的数量）
{
  df<-read.table('03fst/significant/sig.txt',header = T)
  df$pop <- factor(df$pop, levels = c("AFR", "EAS", "EUR", "AMR", "SAS"))
  # 单独
  p<-ggplot(df[df$variable=='total',], aes(x = pop, y = value)) +
    geom_bar(stat = 'identity',fill="#3B4992FF") +
    ylim(0, 120000)+
    labs(x = 'pop', y = 'value', fill = 'variable') +
    theme_bw()+mythem
  ggsave(
    filename = '../graph/sigFST_total.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-ggplot(df[df$variable=='sig',], aes(x = pop, y = value)) +
    geom_bar(stat = 'identity',fill="#BB0021FF") +
    ylim(0, 25000)+
    labs(x = 'pop', y = 'value', fill = 'variable') +
    theme_bw()+mythem
  ggsave(
    filename = '../graph/sigFST_any.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  # 合成一个
  # ggplot(df, aes(x = pop, y = value, fill = variable)) +
  #   geom_bar(stat = 'identity', position = 'dodge') +
  #   ylim(0, 120000)+
  #   labs(x = 'pop', y = 'value', fill = 'variable') +
  #   scale_fill_manual(values = c("#3353A3","#B71C36"))+
  #   theme_bw()+mythem
}
# 富集图
{
  library(tidyr)
  df<-read.table('03fst/significant/mnv_GO.txt',header = T,sep = '\t')
  df<-separate(df,Term,sep = ":",into = c("ID","Term"))
  df<-df[df$Category=='KEGG_PATHWAY' &  as.numeric(df$FDR)<=0.05,]
  p<-ggplot(df,aes(x=Fold.Enrichment,y=Term)) +
    geom_point(aes(size=Count,color=-1*log10(PValue)))+
    scale_colour_gradientn(colors = c("#FFD133", "#FF8D33", "#B0282E"), 
                           values = c(0, 0.5, 1), 
                           guide = "colourbar")+
    labs(
      color='-log10(Pvalue)',
      size="Gene number",
      x="Fold enrichment"
      # y="Pathway name",
      # title="Pathway enrichment")
    )+
    theme_bw()+
    theme(panel.grid=element_blank(),
         plot.margin=unit(rep(2,4),'lines'),
         panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
         title = element_text(vjust=1,size=8,face = "bold", color='black'),
         axis.text.x=element_text(vjust=1,size=12,face = "bold", color='black'),
         axis.title.x=element_text(vjust=0, size=12,face = "bold", color='black'),
         axis.text.y=element_text(vjust=1,size=8,face = "bold", color='black'),
         axis.title.y=element_blank()
         # legend.position = 'none'
    )
  p
  ggsave(
    filename = '../graph/KEGG2.pdf',
    plot = p,width = 6,height = 10,
    units = 'in'
  )

}


