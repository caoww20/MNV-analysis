library(ggplot2)
## Theme
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
# Functions
{
  # Compute Euclidean distance
  euclidean_distance <- function(matrix1, matrix2) {
    sqrt(sum((matrix1 - matrix2)^2))
  }
  # Determine similarity between two matrices
  euclidean_distance_norm <- function(matrix1, matrix2) {
    # First check whether rows and columns are equal
    if (nrow(matrix1)==nrow(matrix2) & ncol(matrix1)==ncol(matrix2)) {
      distance <- euclidean_distance(matrix1, matrix2)
      # Normalize distance into [0,1] by dividing by the max element value.
      # 0 means identical, 1 means completely different; thresholds are adjustable.
      max_value <- max(c(matrix1, matrix2))
      flag <- distance / max_value
      if (flag<=0.1) {
        return(1)
      }else if(flag<=0.5){
        return(2)
      }else{
        return(3)
      }
    }else{
      return(NA)
    }
  }
  # Example
  matrix1 <- matrix(c(1, 2, 3, 4, 6, 7), nrow = 3)
  matrix2 <- matrix(c(2, 3, 4, 5, 6, 7), nrow = 2)
  euclidean_distance_norm(matrix1, matrix2)
  
}
# Overlap of two TF sets
{
  motif_HOCOMOCO<-read.table('./01motif/HOCOMOCO.motif',header = T,stringsAsFactors = F)
  motif_JASPAR<-read.table('./01motif/JASPAR.motif',header = T,stringsAsFactors = F)
  
  x <- list(HOCOMOCO=motif_HOCOMOCO$TF,
            JASPAR=motif_JASPAR$TF)
  
  venn::venn(x,zcolor=2:6,box=T)
}
# Overlap of motifs
{
  
}
# Promoter regions: average number of MNVs (unique MNVs/length after dedup) for all genes, proteins, and lncRNAs
{}
# Distribution of MNV count per gene promoter, e.g., how many genes carry 1, 2, ... MNVs
{}
# Read identified data
{
  # HOCOMOCO TF and motif are one-to-one; JASPAR is not
  readDf<-function(file){
    df<-fread(file ,header=T,stringsAsFactors = F,data.table = F)
    df$TFBSPair<-paste0(df$MNVID,'-',df$motif,'-',df$TF)
    df$TFPair<-paste0(df$MNVID,'-',df$TF)
    df<-df[c(17,18,5,1,2,15,16)]
    df<-df[!duplicated(df,),] # Remove duplicate pairs with different scores (multiple positions)
    return(df)
  }
  hocomoco<-readDf('./02res/hocomoco_f')
  jaspar<-readDf('./02res/jaspar_f')
}
# FIMO results summary: bar plot of gain/loss across two methods
{
  gainLoss<-list()
  gainLoss[['hocomoco']]<-table(hocomoco$Effect)
  gainLoss[['jaspar']]<-table(jaspar$Effect)
  gainLoss<-as.data.frame(do.call(rbind,gainLoss))
  gainLoss$type<-rownames(gainLoss)
  gainLoss<-reshape::melt(gainLoss,id='type')
  
  ggplot(gainLoss, aes(x = type, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Type") +
    ylab("Value") +
    ggtitle("Grouped Bar Chart")+
    theme_bw() + mythem
}

# How many MNVs (xx%) cause changes in TFBS
  print(length(unique(c(hocomoco$MNVID,jaspar$MNVID))))
  
# Total number of unique MNV–TFBS pairs
  length(unique(c(hocomoco$TFBSPair,jaspar$TFBSPair)))
  # 263163

# How many MNV–TFBS (xx%) cause changes different from SNV–TFBS
{
  data<-fread('./02res/mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)  
  table(data$V2)
}
  
# Alluvial plot: compare SNV vs MNV by pairing MNV–TFBS with SNV–TFBS (3x3 category changes)
{
  data<-fread('./02res/snv_mnv_tfbs.txt',header = F,sep = '\t',stringsAsFactors = F,data.table = F)
  # colnames(data)<-c('snv','mnv')
  # table(data) # 3x3 combinations
  # Draw alluvial plot
  library(alluvial)
  data<-as.data.frame(table(data[1:2]),stringsAsFactors=F)
  data<-data[data$Freq>0,]
  colnames(data)<-c('snv','mnv','freq')
  data<-data[order(data$snv),]
  
  set.seed(8) # for nice colors
  cols <- hsv(h = sample(1:8/10), s = sample(3:8)/8, v = sample(3:8)/8)
  ord <- list(NULL, with(data, order(snv, mnv)))
  alluvial(
    data[,1:2], alpha=1,
    freq = data$freq,
    blocks = T,
    gap.width=0.1, # Block spacing
    xw = 0.3,      # Curve bend of lines
    cw = 0.01,     # Label block size
    # col=c('#BDBDBD','#4EAE4A','#984EA3')
    # axis_labels='',
    axes=F,
    ann=F,
    col = cols[match(data$snv, unique(data$mnv))],
    ordering=ord
  )
}

# Upstream 2.5 kb in 100 bp bins
  # CPG curve
  # Curve of different MNV–TFBS counts
  # Boxplot curve of p-values
    # First draw boxplot curve of TFBS p-values for a TF
  # Advanced:
    # Get median p-value per TF in region
    # Then draw boxplot curves for all TFs
  
# demo
  # Use code from paper
  
# MNVQTL enriched in cancer
  