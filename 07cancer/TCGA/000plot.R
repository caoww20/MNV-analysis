library(circlize)
library(stringr)
library(data.table)
# Load themes
{
  mythem<-theme(panel.grid=element_blank(),
                plot.margin=unit(rep(2,4),'lines'),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                title = element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.text.x=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black'),
                legend.position = 'none'
  )
  mythem2<-theme(panel.grid=element_blank(),
                 plot.margin=unit(rep(2,4),'lines'),
                 panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                 title = element_text(vjust=1,size=16,face = "bold", color='black'),
                 axis.text.x=element_text(vjust=0.5,size=16,angle = 45,face = "bold", color='black'),
                 axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                 axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                 axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black'),
                 legend.position = 'none'
  )
}
# Load packages for truncated barplots
{
  library(ggplot2)
  library(plotrix)
  gap.barplot <- function(df, y.cols = 1:ncol(df), sd.cols = NULL, btm = NULL,
                          top = NULL, min.range = 10, max.fold = 5, ratio = 1, gap.width = 1, brk.type = "normal",
                          brk.bg = "white", brk.srt = 135, brk.size = 1, brk.col = "black", brk.lwd = 1,
                          cex.error = 1, ...) {
    if (missing(df))
      stop("No data provided.")
    if (is.numeric(y.cols))
      ycol <- y.cols else ycol <- colnames(df) == y.cols
      if (!is.null(sd.cols))
        if (is.numeric(sd.cols))
          scol <- sd.cols else scol <- colnames(df) == sd.cols
          ## Arrange data
          opts <- options()
          options(warn = -1)
          y <- t(df[, ycol])
          colnames(y) <- NULL
          if (missing(sd.cols))
            sdx <- 0 else sdx <- t(df[, scol])
          sdu <- y + sdx
          sdd <- y - sdx
          ylim <- c(0, max(sdu) * 1.05)
          ## Auto-detect break (btm/top) when not provided
          if (is.null(btm) | is.null(top)) {
            autox <- .auto.breaks(dt = sdu, min.range = min.range, max.fold = max.fold)
            if (autox$flag) {
              btm <- autox$btm
              top <- autox$top
            } else {
              xx <- barplot(y, beside = TRUE, ylim = ylim, ...)
              if (!missing(sd.cols))
                errorbar(xx, y, sdu - y, horiz = FALSE, cex = cex.error)
              box()
              return(invisible(xx))
            }
          }
          ## Set up virtual y limits
          halflen <- btm - ylim[1]
          xlen <- halflen * 0.1 * gap.width
          v_tps1 <- btm + xlen  # virtual top positions
          v_tps2 <- v_tps1 + halflen * ratio
          v_ylim <- c(ylim[1], v_tps2)
          r_tps1 <- top  # real top positions
          r_tps2 <- ylim[2]
          ## Rescale data
          lmx <- summary(lm(c(v_tps1, v_tps2) ~ c(r_tps1, r_tps2)))
          lmx <- lmx$coefficients
          sel1 <- y > top
          sel2 <- y >= btm & y <= top
          y[sel1] <- y[sel1] * lmx[2] + lmx[1]
          y[sel2] <- btm + xlen/2
          sel1 <- sdd > top
          sel2 <- sdd >= btm & sdd <= top
          sdd[sel1] <- sdd[sel1] * lmx[2] + lmx[1]
          sdd[sel2] <- btm + xlen/2
          sel1 <- sdu > top
          sel2 <- sdu >= btm & sdu <= top
          sdu[sel1] <- sdu[sel1] * lmx[2] + lmx[1]
          sdu[sel2] <- btm + xlen/2
          ## bar plot
          xx <- barplot(y, beside = TRUE, ylim = v_ylim, axes = FALSE, names.arg = NULL,
                        ...)
          ## error bars
          if (!missing(sd.cols))
            errorbar(xx, y, sdu - y, horiz = FALSE, cex = cex.error)
          ## Real ticks and labels
          brks1 <- pretty(seq(0, btm, length = 10), n = 4)
          brks1 <- brks1[brks1 >= 0 & brks1 < btm]
          brks2 <- pretty(seq(top, r_tps2, length = 10), n = 4)
          brks2 <- brks2[brks2 > top & brks2 <= r_tps2]
          labx <- c(brks1, brks2)
          ## Virtual ticks
          brks <- c(brks1, brks2 * lmx[2] + lmx[1])
          axis(2, at = brks, labels = labx)
          box()
          ## break marks
          pos <- par("usr")
          xyratio <- (pos[2] - pos[1])/(pos[4] - pos[3])
          xlen <- (pos[2] - pos[1])/50 * brk.size
          px1 <- pos[1] - xlen
          px2 <- pos[1] + xlen
          px3 <- pos[2] - xlen
          px4 <- pos[2] + xlen
          py1 <- btm
          py2 <- v_tps1
          rect(px1, py1, px4, py2, col = brk.bg, xpd = TRUE, border = brk.bg)
          x1 <- c(px1, px1, px3, px3)
          x2 <- c(px2, px2, px4, px4)
          y1 <- c(py1, py2, py1, py2)
          y2 <- c(py1, py2, py1, py2)
          px <- .xy.adjust(x1, x2, y1, y2, xlen, xyratio, angle = brk.srt * pi/90)
          if (brk.type == "zigzag") {
            x1 <- c(x1, px1, px3)
            x2 <- c(x2, px2, px4)
            if (brk.srt > 90) {
              y1 <- c(y1, py2, py2)
              y2 <- c(y2, py1, py1)
            } else {
              y1 <- c(y1, py1, py1)
              y2 <- c(y2, py2, py2)
            }
          }
          if (brk.type == "zigzag") {
            px$x1 <- c(pos[1], px2, px1, pos[2], px4, px3)
            px$x2 <- c(px2, px1, pos[1], px4, px3, pos[2])
            mm <- (v_tps1 - btm)/3
            px$y1 <- rep(c(v_tps1, v_tps1 - mm, v_tps1 - 2 * mm), 2)
            px$y2 <- rep(c(v_tps1 - mm, v_tps1 - 2 * mm, btm), 2)
          }
          par(xpd = TRUE)
          segments(px$x1, px$y1, px$x2, px$y2, lty = 1, col = brk.col, lwd = brk.lwd)
          options(opts)
          par(xpd = FALSE)
          invisible(xx)
  }
  ## Draw error bars
  errorbar <- function(x, y, sd.lwr, sd.upr, horiz = FALSE, cex = 1, ...) {
    if (missing(sd.lwr) & missing(sd.upr))
      return(NULL)
    if (missing(sd.upr))
      sd.upr <- sd.lwr
    if (missing(sd.lwr))
      sd.lwr <- sd.upr
    if (!horiz) {
      arrows(x, y, y1 = y - sd.lwr, length = 0.1 * cex, angle = 90, ...)
      arrows(x, y, y1 = y + sd.upr, length = 0.1 * cex, angle = 90, ...)
    } else {
      arrows(y, x, x1 = y - sd.lwr, length = 0.1 * cex, angle = 90, ...)
      arrows(y, x, x1 = y + sd.upr, length = 0.1 * cex, angle = 90, ...)
    }
  }
  .xy.adjust <- function(x1, x2, y1, y2, xlen, xyratio, angle) {
    xx1 <- x1 - xlen * cos(angle)
    yy1 <- y1 + xlen * sin(angle)/xyratio
    xx2 <- x2 + xlen * cos(angle)
    yy2 <- y2 - xlen * sin(angle)/xyratio
    return(list(x1 = xx1, x2 = xx2, y1 = yy1, y2 = yy2))
  }
  ## Auto-calculate break positions
  .auto.breaks <- function(dt, min.range, max.fold) {
    datax <- sort(as.vector(dt))
    flags <- FALSE
    btm <- top <- NULL
    if (max(datax)/min(datax) < min.range)
      return(list(flag = flags, btm = btm, top = top))
    m <- max(datax)
    btm <- datax[2]
    i <- 3
    while (m/datax[i] > max.fold) {
      btm <- datax[i]
      flags <- TRUE
      i <- i + 1
    }
    if (flags) {
      btm <- btm + 0.05 * btm
      x <- 2
      top <- datax[i] * (x - 1)/x
      while (top < btm) {
        x <- x + 1
        top <- datax[i] * (x - 1)/x
        if (x > 100) {
          flags <- FALSE
          break
        }
      }
    }
    return(list(flag = flags, btm = btm, top = top))
  }
}
# Load base data
{
  # Count MNVs per cancer type
  {
    allMNV<-read.table('./01data/TCGA.mnv',header = F,sep = '\t',stringsAsFactors = F)
    colnames(allMNV)<-c('chr','pos','mnvid','ref','alt','rsid','distance','mnv_type','info')
    cancer_list<-c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC',
                   'KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')
    cancer_num<-list()
    for (i in cancer_list) {
      cancer_num[[i]]=sum(grepl(pattern = i,allMNV$info))
    }
    cancer_num<-do.call(rbind,cancer_num)
  }
  # Count eQTL-MNVs per cancer type
  {
    eQTLMNV<-read.table('./01data/eQTL_fix.txt',header = F,sep = '\t',stringsAsFactors = F)
    colnames(eQTLMNV)<-c('cancer','cis_trans','mnvid','chr','pos','ref','alt','rsid','eQTL_gene','eQTL_symbol','gene_pos','t','beta','st','r','p','fdr','anno','info')
    eQTLMNV<-eQTLMNV[eQTLMNV$fdr<=0.01,]
    tmp<-eQTLMNV[c(1,3)]
    tmp<-tmp[!duplicated(tmp),]
    # 多少个MNV是eQTL-MNV
    length(unique(tmp$mnvid))
    cancer_eQTL_num<-as.data.frame(table(tmp$cancer))
    eQTLMNV_cancer_num<-as.data.frame(table(tmp$mnvid))
    # 10种以上癌症中存在
    nrow(eQTLMNV_cancer_num[eQTLMNV_cancer_num$Freq>=10,])
  }
}
# Overlap with GTEx
{
  1-183129/224158 # 18%
}
# Plot MNV counts and eQTL counts across cancers
{
  df <- data.frame(type = c(rep('mnv',length(cancer_num)),rep('eQTL-mnv',nrow(cancer_eQTL_num))), 
                   cancer = c(rownames(cancer_num),levels(cancer_eQTL_num$Var1)),
                   num=c(cancer_num,cancer_eQTL_num$Freq))
  
  # eQTL counts vary a lot; only plot MNV counts here
  p<-ggplot(df[df$type=='mnv',], aes(cancer,num,fill=type)) +
    # geom_col(position = "dodge") +
    geom_col(position = position_dodge(width = 0.8),fill='#008280FF',color="white") +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(breaks = seq(0, 150000, 50000), limits = c(0, 150000)) +
    theme_bw() + mythem2
  ggsave(
    filename = './graph/cancer_mnv.pdf',
    plot = p,width = 10,height = 4,
    units = 'in'
  )
}
# Plot number of cancers affected per eQTL-MNV
{
  res<-eQTLMNV[c(1,3)]
  res<-res[!duplicated(res),]
  res<-as.data.frame(table(table(res$mnvid)),stringsAsfactor=F)
  res[10,2]<-sum(res[10:nrow(res),2])
  res<-res[1:10,]
  res[1,2]/sum(res[,2]) # 45.72% 唯一组织
  res[10,2]/sum(res[,2]) # 6.61% 泛癌
  # Build color palette
  colors <- colorRampPalette(c("#2A5888", "#C4F5FC"))(10)
  # Draw gradient pie chart
  pie(res$Freq,col = colors)
  # 6*6
}
# Plot circos diagrams https://jokergoo.github.io/circlize_book/book/
{
  mydf<-eQTLMNV[c(3,4,5,18)]
  mydf<-mydf[!duplicated(mydf),]
  mydf<-merge(mydf,eQTLMNV_cancer_num,by.x='mnvid',by.y='Var1',all.x = T,sort =F)
  mydf<-cbind(mydf,str_split_fixed(mydf$pos,',',n=Inf)[,1:2])
  colnames(mydf)[c(6,7)]<-c('start','end')
  # Get promoter
  flag_p<-grep('up:',mydf$anno)
  # Get UTR5
  flag_UTR5<-grep('UTR5',mydf$anno)
  # UTR5/up
  # grep('UTR5|up:',mydf$anno)
  # CDS
  flag_CDS<-grep('exon',mydf$anno)
  # splice
  flag_s<-grep('splice',mydf$anno)
  # UTR3
  flag_UTR3<-grep('UTR3',mydf$anno)
  # intron
  flag_intro<-grep('intro',mydf$anno)
  # other
  flag_other<-setdiff(1:nrow(mydf), unique(c(flag_p,flag_UTR5,flag_CDS,flag_s,flag_UTR3,flag_intro)))
  
  # Plot counts per genomic region
  res <- data.frame(type = c('splice','cds','utr3','utr5','promotor','intro','other'), 
                    num = c(length(flag_s),length(flag_CDS),length(flag_UTR3),length(flag_UTR5),length(flag_p),length(flag_intro),length(flag_other)))
  res$flag<-1:nrow(res)
  gap.barplot(res,y.cols = 'num',btm = 1500,top = 10000,ratio = 1,gap.width = 0.5 , brk.size = 0.4,col='#631879')
  # 6*6
  
  # 使用circ
  {
    # 简化结果
    mydf<-mydf[c(2,6,7,5)]
    mydf$start<-as.integer(mydf$start)
    mydf$end<-as.integer(mydf$end)
    mydf$chr<-paste0('chr',mydf$chr)
    ## 整个 circo
    {
      mycol=c('#9B5400','#435290','#008280FF','#AF282F','#631879FF','#CCCCCC')
  # Outer track: chr labels
      circos.initializeWithIdeogram(species='hg38',chromosome.index = paste0("chr",1:22))
  # Outer track: colored blocks
      set.seed(123)
      circos.initializeWithIdeogram(species='hg19',chromosome.index = paste0("chr",1:22),plotType = NULL)
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        chr = gsub('chr','',CELL_META$sector.index)
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1,luminosity = "bright"))
        circos.text(mean(xlim), mean(ylim)+1.2, chr, cex = 0.7, col = "black",
                    facing = "inside", niceFacing = TRUE)
      }, track.height = 0.06, bg.border = NA)
      
      
  # 'other'
      circos.genomicTrackPlotRegion(mydf[flag_other,], ylim = c(0, 33),
                                    panel.fun = function(region, value, ...) {
                                      circos.genomicPoints(region, value,cex=0.4,pch = 16,col='#CCCCCC')
                                    })
      
      circos.clear()
      
    }
  ## Custom half-circos
    {
      mycol=c('#AF282F','#631879FF','#9B5400','#435290','#008280FF','#ADD3A5','#CCCCCC')
      
  # Duplicate to make two halves symmetric
      human_cytoband = read.cytoband(species = "hg19")$df
      mouse_cytoband = read.cytoband(species = "hg19")$df
      
      mouse_cytoband[ ,1] = paste0("m_", mouse_cytoband[, 1])
      
      cytoband = rbind(human_cytoband, mouse_cytoband)
      head(cytoband)
      chromosome.index = c(rev(paste0("m_chr", c(1:22))),
                           paste0("chr", c(1:22)))
  # Option: directly draw ideograms for both sets
      # circos.initializeWithIdeogram(cytoband, chromosome.index = chromosome.index)
      
  # Only draw chromosome labels
      # circos.par(gap.after = c(rep(1, 21), 5, rep(1, 21), 5))
      circos.initializeWithIdeogram(cytoband, plotType = NULL, 
                                    chromosome.index = chromosome.index)
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
                    gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
      }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
  # Option: draw ideograms only
      # circos.genomicIdeogram(cytoband)
  # Skip ideograms: draw an empty ring and color later in Illustrator
      circos.genomicTrackPlotRegion(mydf, track.height = 0.05)
      # 'other'
      circos.genomicTrackPlotRegion(mydf[unique(flag_other,flag_intro),], ylim = c(0, 33),
                                    panel.fun = function(region, value, ...) {
                                      circos.genomicPoints(region, value,cex=0.3,pch = 16,col='#CCCCCC')
                                    })
  # promoter & UTR5
      circos.genomicTrackPlotRegion(mydf[unique(flag_p,flag_UTR5),],track.index=3,
                                    panel.fun = function(region, value, ...) {
                                      circos.genomicPoints(region, value,cex=0.6,pch = 16,col=alpha('#008280FF',0.5),alpha=0.2)
                                    })
      
  # UTR3
      circos.genomicTrackPlotRegion(mydf[flag_UTR3,],track.index=3,
                                    panel.fun = function(region, value, ...) {
                                      circos.genomicPoints(region, value,cex=0.6,pch = 16,col=alpha('#435290',0.5),alpha=0.2)
                                    })
  # CDS & splice
      circos.genomicTrackPlotRegion(mydf[unique(flag_CDS,flag_s),],track.index=3,
                                    panel.fun = function(region, value, ...) {
                                      circos.genomicPoints(region, value,cex=0.6,pch = 16,col=alpha('#AF282F',0.5),)
                                    })
      circos.clear()
  # A4
    }
  }
  ## Manhattan plots
  {
    library(qqman)
    mydf<-mydf[c(1,2,6,7,5)]
    mydf$start<-as.integer(mydf$start)
    mydf$end<-as.integer(mydf$end)
    colnames(mydf)<-c('MNV','CHR','START','END','Value')
    manhattan(mydf[unique(flag_other,flag_intro),],chr="CHR",bp="START",p="Value",snp="MNV", col=c("#CCCCCC"),chrlabs = as.character(c(1:22)),logp=F,suggestiveline = F, 
              genomewideline = F,ylab="Value",ylim=c(0,33),font.lab=4,cex.lab=1.2,main="",cex=0.8)
    manhattan(mydf[unique(flag_p,flag_UTR5),],chr="CHR",bp="START",p="Value",snp="MNV", col=c("#008280FF"),chrlabs = as.character(c(1:22)),logp=F,suggestiveline = F, 
              genomewideline = F,ylab="Value",ylim=c(0,33),font.lab=4,cex.lab=1.2,main="",cex=0.8)
    manhattan(mydf[flag_UTR3,],chr="CHR",bp="START",p="Value",snp="MNV", col=c("#435290"),chrlabs = as.character(c(1:22)),logp=F,suggestiveline = F, 
              genomewideline = F,ylab="Value",ylim=c(0,33),font.lab=4,cex.lab=1.2,main="",cex=0.8)
    tmp<-rbind(mydf[unique(flag_CDS,flag_s),],c(0,20,1,1,0))
    manhattan(tmp,chr="CHR",bp="START",p="Value",snp="MNV", col=c("#AF282F"),chrlabs = as.character(c(1:22)),logp=F,suggestiveline = F, 
              genomewideline = F,ylab="Value",ylim=c(0,33),font.lab=4,cex.lab=1.2,main="",cex=0.8)
    # 8*3
    
    
    library(CMplot)
    CMplot(mydf[unique(flag_other,flag_intro),c(1,2,3,5)],plot.type="m",LOG10=F,amplify = F,
           threshold=NULL,chr.den.col=NULL,file="pdf",dpi=300,band=3,main= "CMplot intro",file.name='CMplot_intro',col = "#CCCCCC")  
    CMplot(mydf[unique(flag_p,flag_UTR5),c(1,2,3,5)],plot.type="m",LOG10=F,amplify = F,
           threshold=NULL,chr.den.col=NULL,file="pdf",dpi=300,band=3,main= "CMplot utr5",file.name='CMplot_utr5',col = "#008280FF") 
    CMplot(mydf[flag_UTR3,c(1,2,3,5)],plot.type="m",LOG10=F,amplify = F,
           threshold=NULL,chr.den.col=NULL,file="pdf",dpi=300,band=3,main= "CMplot utr3",file.name='CMplot_utr3',col = "#435290") 
    tmp<-rbind(mydf[unique(flag_CDS,flag_s),],c(0,20,1,1,0))
    CMplot(tmp[,c(1,2,3,5)],plot.type="m",LOG10=F,amplify = F,
           threshold=NULL,chr.den.col=NULL,file="pdf",dpi=300,band=3,main= "CMplot cds",file.name='CMplot_cds',col = "#AF282F") 

     
  }
}
# Plot beta values across regions
{
  # Filter extreme values
  filterExtremum<-function(df){
    colnames(df)<-c('type','num')
    x=df$num
    iqr <- IQR(x)
    upper_bound <- median(x) + 1.5 * iqr
    lower_bound <- median(x) - 1.5 * iqr
    df_filtered <- df[df$num >= lower_bound & df$num <= upper_bound,]
    return(df_filtered)
  }
  
  # Use eQTL-MNV dataset
  # Get promoter
  flag_p<-grep('up:',eQTLMNV$anno)
  # Get UTR5
  flag_UTR5<-grep('UTR5',eQTLMNV$anno)
  # CDS
  flag_CDS<-grep('exon',eQTLMNV$anno)
  # splice
  flag_s<-grep('splice',eQTLMNV$anno)
  # UTR3
  flag_UTR3<-grep('UTR3',eQTLMNV$anno)
  # intron
  flag_intro<-grep('intro',eQTLMNV$anno)
  # other
  flag_other<-setdiff(1:nrow(eQTLMNV), unique(c(flag_p,flag_UTR5,flag_CDS,flag_s,flag_UTR3,flag_intro)))
  
  res <- data.frame(type = c(rep('splice',length(flag_s)),rep('cds',length(flag_CDS)),rep('utr3',length(flag_UTR3)),rep('utr5',length(flag_UTR5)),rep('promotor',length(flag_p)),rep('intro',length(flag_intro)),rep('other',length(flag_other))), 
                    num = c(eQTLMNV[flag_s,'beta'],eQTLMNV[flag_CDS,'beta'],eQTLMNV[flag_UTR3,'beta'],eQTLMNV[flag_UTR5,'beta'],eQTLMNV[flag_p,'beta'],eQTLMNV[flag_intro,'beta'],eQTLMNV[flag_other,'beta']))
  res$num<-abs(res$num)
  res$type <- factor(res$type,levels=c("splice","cds","utr5","promotor","utr3","intro",'other'))
  
  ggplot(res,aes(x=type,y=num))+
    geom_boxplot(aes(fill = type),outlier.shape = NA)+
    ylim(0, 0.75)+
    labs(x='',y='beta')+theme_bw()+mythem
  # Exclude splice in violin plot due to low counts
  library(smplot2)
  p<-ggplot(res[res$type!='splice',], aes(x = type, y = num, fill = type))+
    geom_half_violin(aes(fill = type),
                     position = position_nudge(x = .15, y = 0),
                     adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
    geom_boxplot(aes(x = type,y = num, fill = type),
                 outlier.shape = NA,
                 width = .1,
                 color = "black")+
    labs(x=NULL,y='The Beta of eQTL-MNVs')+
    scale_fill_manual(values = colorRampPalette(c("#2A5888", "#C4F5FC"))(6))+
    scale_y_continuous(breaks=seq(0, 0.75, by = 0.25),limits = c(0, 0.75))+
    theme_bw()+mythem
  ggsave(
    filename = './graph/beta.pdf',
    plot = p,width = 6,height = 5,
    units = 'in'
  )
  # Alternative: filter outliers before plotting (discarded)
  tmp=list()
  for (i in c("splice","cds","utr3","utr5","promotor","intro",'other')) {
    tmp[[i]]<-filterExtremum(res[res$type==i,])
  }
  res<-do.call(rbind,tmp)
  res$type <- factor(res$type,levels=c("splice","cds","utr3","utr5","promotor","intro",'other'))
  ggplot(res,aes(x=type,y=num))+
    geom_boxplot(aes(fill = type),outlier.shape = NA)+
    ylim(0, 0.75)+
    labs(x='',y='beta')+theme_bw()+mythem
}
# Overlap with cancer genes
{
  library(eulerr)
  # Cancer genes
  cancer_gene1<-read.table('/data/jinww/04reference/publicDB/cancer_gene/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv',header = T,sep = '\t',stringsAsFactors = F)
  cancer_gene2<-read.table('/data/jinww/04reference/publicDB/cancer_gene/NCG_cancerdrivers_annotation_supporting_evidence.tsv',header = T,sep = '\t',stringsAsFactors = F)
  cancer_gene<-unique(c(cancer_gene1$SYMBOL,cancer_gene2$symbol))
  # Overlap between eQTL-MNV-associated genes and cancer genes
  {
    eQTL_gene<-unique(eQTLMNV$eQTL_symbol)
    over_gene<-intersect(eQTL_gene, cancer_gene)
    
    fit1 <- euler(c("A" = length(cancer_gene)-length(over_gene), "B" = length(eQTL_gene)-length(over_gene), 
                    "A&B" = length(over_gene)))
    
    fill1=list(fill=c("#2A5888","#89CEED"),ncol=2,
               alpha=1)
    
    plot(fit1,fills = fill1,
         #labels = list(col="white",font=3,cex=2),
         quantities = list(col="black",cex=2),
         legend = list(labels=c("cancer","mnv"),cex=2,col="white",fontsize=10,side='top',pch=21),
         edges = list(lwd=3)
    )
  # 6*6
    tmp<-eQTLMNV[eQTLMNV$eQTL_symbol%in%over_gene,]
    nrow(tmp)
  }
  # Count eQTL-MNVs located on cancer genes
  {
    mnv_gene<-read.table('01data/eQTL_fix.gene',header = F,stringsAsFactors = F)
    mnv_gene<-mnv_gene$V1
    over_gene<-intersect(cancer_gene, mnv_gene)
    
    fit1 <- euler(c("A" = length(cancer_gene)-length(over_gene), "B" = length(mnv_gene)-length(over_gene), 
                    "A&B" = length(over_gene)))
    
    fill1=list(fill=c("#e6194B","#5B9BD5"),ncol=2,
               alpha=1)
    
    plot(fit1,fills = fill1,
         #labels = list(col="white",font=3,cex=2),
         quantities = list(col="black",cex=2),
         legend = list(labels=c("cancer","mnv"),cex=2,col="white",fontsize=10,side='top',pch=21),
         edges = list(lwd=3,lty=2)
    )
  }
  
}
# LUAD results
{
  eQTLGwas<-read.table('01data/eQTL_gwas.txt',header = F,stringsAsFactors = F,sep = '\t')
  # eQTLGwas<-na.omit(eQTLGwas)
  rownames(eQTLGwas)<-NULL
  colnames(eQTLGwas)<-c('cancer','cis_trans','mnvid','chr','pos','ref','alt','rsid','eQTL_gene',
                        'eQTL_symbol','gene_pos','t','beta','st','r','p','fdr','gwas_beta','gwas_p','anno','info')
  # lung_cancer_gene
  lung_cancer_gene<-unique(c(cancer_gene1[cancer_gene1$CANCER_TYPE=='LUAD' |cancer_gene1$CANCER_TYPE=='LUSC','SYMBOL'],cancer_gene2[grepl('lung',cancer_gene2$primary_site),'symbol']))
  # eQTL gene is a lung cancer gene
  a=which(eQTLGwas$eQTL_symbol%in%lung_cancer_gene)
  # MNV annotation contains lung cancer gene name
  b=c()
  for (i in intersect(lung_cancer_gene, mnv_gene)) {
    if (i=='NF1' | i=='APC') {
      next
    }
    b=c(b,grep(i,eQTLGwas$anno))
  }
  b=unique(b)
  
  a_b=intersect(a,b)
  ab=unique(c(a,b))
  res<-eQTLGwas[ab,]
  
}