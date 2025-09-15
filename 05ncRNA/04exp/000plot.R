# Themes
{
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(tidyr)
  mythem=theme(panel.grid=element_blank(),
               plot.margin=unit(rep(2,4),'lines'),
               # legend.position="none",
               panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
               # linewidth replaces deprecated element_rect
               # title = element_text(vjust=1,size=16,face = "bold", color='black'),
               axis.text.x=element_text(vjust=1,size=8,face = "bold", colour='black'),
               axis.title.x=element_text(vjust=0, size=8,face = "bold", color='black'),
               axis.text.y=element_text(vjust=0.5,size=8,face = "bold", color = 'black'),
               axis.title.y=element_text(vjust=2, size=8,face = "bold", color='black')
  )
}
# Helpers
{
  # Find closest upper multiple
  find_closest_multiple <- function(vector) {
    max_value <- max(vector)
    remainder <- max_value %% 50
    
    if (remainder == 0) {
      closest_multiple <- max_value
    } else {
      closest_multiple <- max_value + (50 - remainder)
    }
    
    return(closest_multiple)
  }
  # Find closest upper multiple (smarter bucketing)
  find_closest_multiple <- function(vector) {
    max_value <- max(vector)
    
    if (max_value < 20) {
      remainder <- max_value %% 5
      closest_multiple <- max_value + (5 - remainder)
    } else if (max_value < 50) {
      remainder <- max_value %% 10
      closest_multiple <- max_value + (10 - remainder)
    } else {
      remainder <- max_value %% 50
      if (remainder == 0) {
        closest_multiple <- max_value
      } else {
        closest_multiple <- max_value + (50 - remainder)
      }
    }
    
    return(closest_multiple)
  }
  # t.test comparison
  caltt<-function(x){
    # print(x[1])
    # p_value <- wilcox.test(x[1:3], x[4:6])$p.value
    if(x[1]==x[2] & x[1]==x[3] & x[4]==x[5] & x[4]==x[6]){
      return(1)
    }
    # Avoid: Error in t.test.default(...): data are essentially constant
    p_value <- t.test(x[1:3], x[4:6])$p.value
    return(p_value)
  }
  # Fold change (log2)
  calDiff<-function(x){
    diff<-log(mean(x[1:3]+0.01)/mean(x[4:6]+0.01),2)
    return(diff)
  }
  # Differential expression filter
  calSig<-function(df){
    flag<-rowSums(df)
    df<-df[which(flag>5),]
    df$p<-apply(df, 1, caltt)
    # df$fdr<-p.adjust(df$p) # With few samples, multiple testing may be unnecessary
    df$foldchange2<-apply(df, 1, calDiff)
    df<-df[abs(df$foldchange2)>log2(1.2) & df$p<0.05,]
    return(df)
  }
}
# Base data
{
  gene<-read.table('02res/symbol.map',header = F,stringsAsFactors = F,sep = '\t')
  colnames(gene)<-c('trans','gene','symbol')
  sample<-read.table('02res/sample.info',header = T,stringsAsFactors = F,sep = '\t')
}
# Gain/Loss counts for 7 miRNAs affecting seed regions
{
  df<-read.table('../03miRBS_UTR3/05altMIR/merge/UTR3/diff_MNV_vs_ref.txt',header = F,sep = '\t',stringsAsFactors = F)
  
  data<-df[df$V16=='gain',1:2]
  data$V1<-str_split_fixed(data$V1,'[|]',n=2)[,1]
  data<-data[!duplicated(data),]
  res<-merge(data,gene,by.x='V1',by.y='trans',all.x = T,sort = F)
  res<-res[2:3]
  res<-res[!duplicated(res),] 
  gain<-table(res$V2) 
  
  data<-df[df$V16=='loss',1:2]
  data$V1<-str_split_fixed(data$V1,'[|]',n=2)[,1]
  data<-data[!duplicated(data),]
  res<-merge(data,gene,by.x='V1',by.y='trans',all.x = T,sort = F)
  res<-res[2:3]
  res<-res[!duplicated(res),] 
  loss<-table(res$V2)
  
  res<-cbind(gain,loss)
  res<-reshape2::melt(res,by=0)
  res$Var1<-factor(res$Var1,levels = c('hsa-miR-4440','hsa-miR-6826-5p','hsa-miR-564','hsa-miR-4472','hsa-miR-3939','hsa-miR-548l','hsa-miR-302c-5p','hsa-miR-4293','hsa-miR-585-3p','hsa-miR-6726-3p','hsa-miR-449c-3p','hsa-miR-3126-5p','hsa-miR-6787-5p','hsa-miR-593-5p','hsa-miR-6796-3p','hsa-miR-4781-3p','hsa-miR-4514'))
  p<-ggplot(res, aes(x = Var1, y = value, fill = Var2)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("#2A5888","#4396B1")) +
    # ylim(0,3500)+
    # scale_y_continuous(breaks = seq(0, 3500, 500), limits = c(0, 3500), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 3500, 500), limits = c(0, 3500)) +
    labs(x = '', y = "Number of Genes Affected by miRNA") +
    theme_bw()+mythem
  p
  ggsave(
    filename = './graph/seven_miRNA.pdf',
    plot = p,width = 8,height = 4,
# Clustering plots
  ) 
  # Read TPM
}
# Clustering plots
{
  # Read TPM
  # Sample clustering
  # Remove zero-variance rows (uninformative and harmful for PCA)
  data<-data[4:ncol(data)]
  
  # Optional preprocessing
  # MAD (median absolute deviation) to measure expression variability
  # Consider filtering by top-k MAD genes for downstream analyses
  # 计算中值绝对偏差 (MAD, median absolute deviation)度量基因表达变化幅度
  # 在基因表达中，尽管某些基因很小的变化会导致重要的生物学意义，
  # Choose top-N (here 5000)
  # 可以选择根据mad值过滤，取top 50， 500等做分析
  mads <- apply(data, 1, mad)
  data <- data[rev(order(mads)),]
  # PCA input: rows=samples, cols=variables
  dim(data)
  
  # Pay attention to the format of PCA input
  # Rows are samples and columns are variables
  data_t <- t(data)
    # Read sample metadata
    # sample<-read.table('sample.info',header = F,stringsAsFactors = F,sep = '\t')
    # colnames(sample)<-c('id','type')
  if(sampleFile != "") {
    # Read sample metadata
    # sample<-read.table('sample.info',header = F,stringsAsFactors = F,sep = '\t')
    # colnames(sample)<-c('id','type')
    sample <- read.table(sampleFile,header = T, row.names=1,sep="\t")
    data_t_m <- merge(data_t, sample[c(2)], by=0)
  # PCA
    data_t <- data_t_m[,-1]
  }
  
  # Run PCA on selected variables
  pca <- prcomp(data_t[,1:variableL], scale=T)
  print(str(pca))
  library(factoextra)
  # Exclude problematic sample 2T0008
  fviz_pca_ind(pca, repel=T)
  fviz_pca_ind(pca, col.ind=data_t$type2, mean.point=F, addEllipses = T, legend.title="Groups")
  
  # Exclude problematic sample 2T0008
  data_t_2<-data_t[c(1:7,9:15),]
# Batch correction with ComBat
  fviz_pca_ind(pca, col.ind=data_t_2$type2, mean.point=F, addEllipses = T, legend.title="Groups")
  
}
# Batch correction with ComBat
{
  library(sva)
  # gene
  edata<-read.table('./02res/gene_tpm.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
  rownames(edata)<-edata$gene
  edata<-edata[4:ncol(edata)]
  edata<-edata[apply(edata, 1, mean) >1,]
  
  pheno<-sample[2:3]
  rownames(pheno)<-sample$sample
  pheno$batch<-1
  pheno[8,3]<-2
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE) # parametric
  batch = pheno$batch
  # combat_edata2 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=FALSE, mean.only=TRUE) # non-parametric
  
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE) # parametric
  combat_edata[combat_edata < 0] <- 0
  # combat_edata2 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=FALSE, mean.only=TRUE) # non-parametric
  data_t<-t(combat_edata)
  pca <- prcomp(data_t, scale=T)
  fviz_pca_ind(pca, col.ind=pheno$type2, mean.point=F, addEllipses = T, legend.title="Groups")
  write.table(combat_edata,'./graph/gene_tpm_fix.txt',quote = F,sep = '\t')
  
  #trans
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE) # parametric
  rownames(edata)<-edata$trans
  edata<-edata[4:ncol(edata)]
  edata<-edata[apply(edata, 1, mean) >1,]
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE) # parametric
  combat_edata[combat_edata < 0] <- 0
 
  data_t<-t(combat_edata)
# miRNA expression barplots
  fviz_pca_ind(pca, col.ind=pheno$type2, mean.point=F, addEllipses = T, legend.title="Groups")
  write.table(combat_edata,'./graph/trans_tpm_fix.txt',quote = F,sep = '\t')
}
# miRNA expression plots
{
  miRNA_exp<-read.table('02res/miRNA_RPM.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
  miRNA_exp<-as.data.frame(t(miRNA_exp[1:4,]))
  colnames(miRNA_exp)<-c('WT','SNV1','SNV2','MNV')
  miRNA_exp$treat<-c(rep('WT',3),rep('MNV',3),rep('SNV2',3),rep('SNV1',3),rep('NC',3))
  miRNA_exp<-reshape2::melt(miRNA_exp,id='treat')
  # Plot
  miRNA_exp$value<-log2(miRNA_exp$value+1)
  miRNA_exp$treat <- factor(miRNA_exp$treat,levels=c("NC","WT","MNV","SNV1","SNV2"))
  
  # Plot figure
  p<-ggplot(miRNA_exp, aes(x = treat, y = value,fill=variable)) +
    geom_bar(stat = 'identity', position = "stack", width=0.8) +
    scale_fill_manual(values = c("#DBDBDB","#4396B1","#89CEED","#2A5888")) +
    scale_y_continuous(breaks = seq(0, 15, 5), limits = c(0, 15)) +
    labs(x = 'treat', y = "log2(RPM)") +
    theme_bw()+mythem
  ggsave(
    filename = './graph/miRNA_exp.pdf',
  # Grouped bars (alternate)
    units = 'in'
  )
  
  # Grouped bar chart (less preferred visually)
  ggplot(miRNA_exp, aes(x = treat, y = value, fill = variable)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("#DBDBDB","#4396B1","#89CEED","#2A5888")) +
# KEGG
    labs(x = 'treat', y = "log2(RPM)") +
    theme_bw()+mythem
}
# KEGG
{
  # MNV_vs_WT
  # gene
  {
    # data<-read.table('./02res/gene_tpm.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
    # rownames(data)<-data$gene
    # data<-data[4:ncol(data)]
    data<-read.table('./graph/gene_tpm_fix.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
    MNV_vs_WT<-calSig(data[c(4:6,1:3)])
  # MNV_vs_NC<-calSig(data[c(4:6,13:15)])
  # WT_vs_NC<-calSig(data[c(1:3,13:15)])
    write.table(unique(gene[which(gene$gene%in%rownames(MNV_vs_WT)),3]),'./graph/MNV_vs_WT_gene.diff',quote = F,row.names = F,col.names = F)
    # write.table(unique(gene[which(gene$gene%in%rownames(MNV_vs_NC)),3]),'MNV_vs_NC.diff',quote = F,row.names = F,col.names = F)
    # write.table(unique(gene[which(gene$gene%in%rownames(WT_vs_NC)),3]),'WT_vs_NC.diff',quote = F,row.names = F,col.names = F)
  }
  # trans
  {
    # data<-read.table('./02res/trans_tpm.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
    # rownames(data)<-data$trans
    # data<-data[4:ncol(data)]
    data<-read.table('./graph/trans_tpm_fix.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
    MNV_vs_WT<-calSig(data[c(4:6,1:3)])
    write.table(unique(gene[which(gene$trans%in%rownames(MNV_vs_WT)),3]),'./graph/MNV_vs_WT_trans.diff',quote = F,row.names = F,col.names = F)
  }
  # KEGG enrichment MNV_vs_WT
  {
    # kegg<-read.table('./graph/MNV_vs_WT_gene_fix.txt',sep = '\t',quote = '',header = TRUE)
    kegg<-read.table('./graph/MNV_vs_WT_trans_fix.txt',sep = '\t',quote = '',header = TRUE)
    kegg<-kegg[kegg$Category=='KEGG_PATHWAY' & kegg$FDR<=0.05,]
    kegg<-separate(kegg,Term,sep = ":",into = c("ID","Term"))

    p<-ggplot(kegg,aes(x=Fold.Enrichment,y=ID)) +
      geom_point(aes(size=Count,color=-1*log10(PValue)))+
      scale_colour_gradientn(colors = c("#FFD133", "#FF8D33", "#B0282E"))+
      labs(
        color='-Log10(Pvalue)',
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
      filename = './graph/KEGG_trans_ID.pdf',
      plot = p,width = 6,height = 10,
      units = 'in'
    )
  }
}
# Differential analysis
{
  # gene
  {
    data<-read.table('./graph/gene_tpm_fix.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE,row.names = 1)
    # rownames(data)<-data$gene
    # data<-data[4:ncol(data)]
    # data$`02T0008`<-rowMeans(data[c(7,9)])
    data_gene<-data
    MNV_vs_WT<-calSig(data[c(4:6,1:3)])
    MNV_vs_NC<-calSig(data[c(4:6,13:15)])
    MNV_vs_SNV1<-calSig(data[c(4:6,7:9)])
    MNV_vs_SNV2<-calSig(data[c(4:6,10:12)])
    venn::venn(list(A=rownames(MNV_vs_WT),B=rownames(MNV_vs_NC),C=rownames(MNV_vs_SNV1),D=rownames(MNV_vs_SNV2)),zcolor= 2:5,box=F,sncs=2,ilcs=1.5)
    ensembl_id <- Reduce(intersect, list(A=rownames(MNV_vs_WT),B=rownames(MNV_vs_NC),C=rownames(MNV_vs_SNV1),D=rownames(MNV_vs_SNV2)))
    common_elements<-unique(gene[which(gene$gene%in%ensembl_id),3])
    write.table(common_elements,'./graph/common_diff.gene',quote = F,row.names = F,col.names = F)
  # 600*600
  }
  # trans
  {
    data<-read.table('./graph/trans_tpm_fix.txt',header = T,stringsAsFactors = F,sep = '\t', check.names = FALSE)
    # rownames(data)<-data$trans
    # data<-data[4:ncol(data)]
    # data$`02T0008`<-rowMeans(data[c(7,9)])
    data_trans<-data
    # data_trans <- data_trans[apply(data_trans, 1, var)!=0,]
    MNV_vs_WT<-calSig(data[c(4:6,1:3)])
    MNV_vs_NC<-calSig(data[c(4:6,13:15)])
    MNV_vs_SNV1<-calSig(data[c(4:6,7:9)])
    MNV_vs_SNV2<-calSig(data[c(4:6,10:12)])
    # trans
    venn::venn(list(A=rownames(MNV_vs_WT),B=rownames(MNV_vs_NC),C=rownames(MNV_vs_SNV1),D=rownames(MNV_vs_SNV2)),zcolor= 2:5,box=F,sncs=2,ilcs=1.5)
    common_elements_trans <- Reduce(intersect, list(A=rownames(MNV_vs_WT),B=rownames(MNV_vs_NC),C=rownames(MNV_vs_SNV1),D=rownames(MNV_vs_SNV2)))
    write.table(common_elements_trans,'./graph/common_diff.trans',quote = F,row.names = F,col.names = F)
    # map to gene
    venn::venn(list(A=gene[which(gene$trans%in%rownames(MNV_vs_WT)),3],B=gene[which(gene$trans%in%rownames(MNV_vs_NC)),3],
                    C=gene[which(gene$trans%in%rownames(MNV_vs_SNV1)),3],D=gene[which(gene$trans%in%rownames(MNV_vs_SNV2)),3]),zcolor= 2:5,box=F,sncs=2,ilcs=1.5)
    # Note: Venn diagram needs manual adjustment if a symbol is '.'
    common_elements <- Reduce(intersect, list(A=gene[which(gene$trans%in%rownames(MNV_vs_WT)),3],B=gene[which(gene$trans%in%rownames(MNV_vs_NC)),3],
                                              C=gene[which(gene$trans%in%rownames(MNV_vs_SNV1)),3],D=gene[which(gene$trans%in%rownames(MNV_vs_SNV2)),3]))
    write.table(common_elements,'./graph/common_diff.trans_to_gene',quote = F,row.names = F,col.names = F)
    # 600*600
  }
  # Note: trans may raise an error
  # Error in t.test.default(x[1:3], x[4:6]) : data are essentially constant
  # Due to identical values in a group, e.g., t.test(c(0.1,0.1,0.1),c(0.2,0.2,0.2))
}
# Get transcripts with differential binding
{
  diff_bind<-read.table('./02res/diff.targetscan',header = F,stringsAsFactors = F,sep = '\t')
  colnames(diff_bind)<-c('trans','miRNA','start','end','seed','type')
  diff_bind$trans<-str_split_fixed(diff_bind$trans,'[|]',2)[,1]
  diff_bind<-diff_bind[diff_bind$trans%in%common_elements_trans,]
  diff_bind<-merge(diff_bind,gene,by='trans',all.x=T,sort=F)
  write.table(diff_bind,'./graph/sig_tar.res',sep = '\t',quote = F,row.names = F)
}
{
  diff_bind<-read.table('./02res/diff.miranda',header = F,stringsAsFactors = F,sep = '\t')
  diff_bind<-diff_bind[c(1,2,3,4,7,8,12,13,14,15)]
  colnames(diff_bind)<-c('trans','miRNA','score','energy','start','end','fig1','fig2','fig3','type')
  diff_bind$trans<-str_split_fixed(diff_bind$trans,'[|]',2)[,1]
  diff_bind<-diff_bind[diff_bind$trans%in%common_elements_trans,]
  diff_bind<-merge(diff_bind,gene,by='trans',all.x=T,sort=F)
  write.table(diff_bind,'./graph/sig_miranda.res',sep = '\t',quote = F,row.names = F)
}
# Example expression plots
{
  # trans
  for (trans_id in diff_bind$trans) {
    # trans_id='ENST00000320895'
    gene_id=unique(diff_bind[diff_bind$trans==trans_id,'gene'])
    gene_symbol=unique(diff_bind[diff_bind$trans==trans_id,'symbol'])
    
    trans_exp=as.data.frame(cbind(t(data_trans[which(rownames(data_trans)==trans_id),]),sample$type2))
    colnames(trans_exp)<-c('tpm','treat')
    trans_exp$tpm<-as.numeric(trans_exp$tpm)
    trans_exp$treat <- factor(trans_exp$treat,levels=c("NC","WT","MNV","SNV1","SNV2"))
    # Bar colors
    # bar_colors <- c("#1A325F","#2272B4","#9AC9E0","#F5AC89","#BD2A34")
    bar_colors <- c("#DBDBDB","#A0A0A0","#2A5888","#4396B1","#89CEED")
    # Plot
    max_value<-find_closest_multiple(trans_exp$tpm)
    p<-ggplot(trans_exp, aes(x = treat, y = tpm)) +
      geom_bar(stat = "summary", fun = "mean", fill = bar_colors, color = "white",alpha=0.8) + # mean bars
      geom_point(position=position_jitter(w=0.1,h=0),size = 3,shape=1,alpha=0.8) + # jittered points
      scale_y_continuous(breaks = seq(0, max_value, max_value/5), limits = c(0, max_value)) +
      labs(x = paste0(gene_symbol,'_',trans_id), y = "TPM") +
      theme_bw()+mythem+
      scale_fill_manual(values = bar_colors)
    ggsave(
      filename = paste0('./graph/graph_demo_trans/t_',gene_symbol,'_',trans_id,'.pdf'),
      plot = p,width = 6,height = 4,
      units = 'in'
    ) 
  }
  # gene
  for (gene_id in unique(diff_bind$gene)) {
    # gene_id='ENSG00000049656'
    gene_symbol=unique(gene[gene$gene==gene_id,'symbol'])
    gene_exp=as.data.frame(cbind(t(data_gene[which(rownames(data_gene)==gene_id),]),sample$type2))
    # if (ncol(gene_exp)==1) {
    #   next()
    # }
    colnames(gene_exp)<-c('tpm','treat')
    gene_exp$tpm<-as.numeric(gene_exp$tpm)
    gene_exp$treat <- factor(gene_exp$treat,levels=c("NC","WT","MNV","SNV1","SNV2"))
    # Bar colors
    bar_colors <- c("#DBDBDB","#A0A0A0","#2A5888","#4396B1","#89CEED")
    # Plot
    max_value<-find_closest_multiple(gene_exp$tpm)
    p<-ggplot(gene_exp, aes(x = treat, y = tpm)) +
      geom_bar(stat = "summary", fun = "mean", fill = bar_colors, color = "white",alpha=0.8) + # mean bars
      geom_point(position=position_jitter(w=0.1,h=0),size = 3,shape=1,alpha=0.8) + # jittered points
      scale_y_continuous(breaks = seq(0, max_value, max_value/5), limits = c(0, max_value)) +
      labs(x = gene_symbol, y = "TPM") +
      theme_bw()+mythem+
      scale_fill_manual(values = bar_colors)
    ggsave(
      filename = paste0('./graph/graph_demo_gene/',gene_symbol,'.pdf'),
      plot = p,width = 6,height = 4,
      units = 'in'
    ) 
  }
}



