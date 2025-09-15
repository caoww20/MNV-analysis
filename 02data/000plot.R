library(ggplot2)
library(ggsci)
library(scales)
library(data.table)
library(tidyverse)
# Extract AAAS palette
color1 = pal_aaas()(10)
show_col(color1)
colorPalette<-c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF")
colorPalette<-c("#e41a1c","#377eb8","#4daf4a","#9ecae1","#6baed6","#4292c6")
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
# Plotting helpers
{
  graph_venn<-function(mylist,mycolors=c("#BB0021FF","#3B4992FF"),my_alpha=0.5){
    A=unique(mylist[[1]])
    B=unique(mylist[[2]])
    ab=intersect(A,B)
    
    print(length(A))
    print(length(B))
    print(length(ab))
    
    
    # myvector<-c("A" = 707, "B" = 401, 
    #                 "A&B" = 347)
    # mycolors<-c("#e6194B","#5B9BD5")
    # mycolor_n=2
    # my_alpha=0.5
    # mylabels=c("JASPAR","HOCOMOCO")
    
    library(eulerr)
    myvector<-c("A" = length(A)-length(ab), "B" = length(B)-length(ab), "A&B" = length(ab))
    fit1 <- euler(myvector)
    
    fill1=list(fill=mycolors,ncol=2,alpha=my_alpha)
    
    p<-plot(fit1,fills = fill1,
            #labels = list(col="white",font=3,cex=2),
            quantities = list(col="black",cex=2),
            legend = list(labels=c('A','B'),cex=2,col="white",fontsize=10,side='top',pch=21),
            edges = list(lwd=3)
            # edges = list(lwd=3,lty=2)
    )
         
    return(p)
  }
}
# Truncated (gap) barplot utilities
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
          ## Auto-detect break (btm/top) if not provided
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
  ## Auto-calc break positions
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
# Load base dataset
{
  mnv<-fread('02datasets/flag_datasets', sep = "\t", quote = '', header = T, stringsAsFactors = F, na.strings = "", data.table = F)
  # colnames(mnv)<-c('chr','pos','mnvid','ref','alt','rsid','distance','mnvtype','1000G','GTEx','UKB20w','UKB50w','TCGA','GnomAD_WGS','GnomAD_WES','numberDB')
}
# Count by dataset and MNV joint type
{
  mnvdb=list()
  mnvdb[['1000G']]=table(mnv[mnv$`1000G`==1,'mnvtype'])
  mnvdb[['GTEx']]=table(mnv[mnv$GTEx==1,'mnvtype'])
  mnvdb[['UKB20w']]=table(mnv[mnv$UKB20w==1,'mnvtype'])
  mnvdb[['UKB50w']]=table(mnv[mnv$UKB50w==1,'mnvtype'])
  mnvdb[['TCGA']]=table(mnv[mnv$TCGA==1,'mnvtype'])
  mnvdb[['GnomAD_WGS']]=table(mnv[mnv$GnomAD_WGS==1,'mnvtype'])
  mnvdb[['GnomAD_WGS']]=c(mnvdb[['GnomAD_WGS']], rep(0, 4 - length(mnvdb[['GnomAD_WGS']])))
  mnvdb[['GnomAD_WES']]=table(mnv[mnv$GnomAD_WES==1,'mnvtype'])
  mnvdb[['GnomAD_WES']]=c(mnvdb[['GnomAD_WES']], rep(0, 4 - length(mnvdb[['GnomAD_WES']])))
  mnvdb[['total']]=table(mnv$mnvtype)
  mnvdb<-do.call(rbind,mnvdb)
  mnvdb<-cbind(mnvdb,rowSums(mnvdb))
  colnames(mnvdb)<-paste0(1:5,'joint')
  write.table(mnvdb,'graph/mnvdb.txt',sep = '\t',quote = F,row.names = T,col.names = T)
}
# Overlap across datasets (excluding gnomAD)
{
  # mnv<-fread('/data/jinww/mnv/analyse/part0/humanMNV/flag_datasets', sep = "\t", quote = '', header = T, stringsAsFactors = F, na.strings = "", data.table = F)
  # Prepare dataset
  {
    data<-mnv[9:13]
    data<-data[!(rowSums(data)==0),]
    data[which(data$`1000G`==1),1]<-'1000G'
    data[which(data$`GTEx`==1),2]<-'GTEx'
    data[which(data$`UKB20w`==1),3]<-'UKB20w'
    data[which(data$`UKB50w`==1),4]<-'UKB50w'
    data[which(data$`TCGA`==1),5]<-'TCGA'
    data[data==0]<-NA
    getMerge<-function(x){
      x<-na.omit(x)
      y=paste(x,collapse = '&')
      return(y)
    }
    res<-apply(data, 1, getMerge)
    res2<-table(res)
    write.table(res2,'./graph/upsetR.txt',quote = F,row.names = F,col.names = T)
    remove(data)
  }
  # Re-process dataset for plotting
  {
    res<-read.table('./graph/upsetR.txt',header = T,stringsAsFactors = F,sep = ' ')
    res<-res[order(-res$Freq),]
    mydata<-res$Freq
    names(mydata)<-res$res
  # UpSet plot
    library(UpSetR)
    data <- fromExpression(mydata)
    upset(data,
          # sets = colnames(data),# sets to display
          nset = 7, # Limit number of sets
          nintersects = 30, # Number of intersections to plot
          # order.by = c('degree','freq'), decreasing = c(F, T),# ordering
          order.by = c("freq"), decreasing = c(TRUE),
          mb.ratio = c(0.7, 0.3),# Ratio between main bar and matrix
          number.angles = 0,# Angle for numbers on bars
          point.size = 6,# Circle size in matrix
          line.size = 1.5, # Size of lines/points in matrix
          show.numbers = F, # Show numbers on bars
          # shade.color = "red",# Shadow color in matrix
          mainbar.y.label = "Intersection size", # Y-axis label of main barplot
          sets.x.label = "Set Size", # X-axis label of set size bars
          main.bar.color = "#631879FF", # Color of main bars
          sets.bar.color = "#008280FF",# Color of set size bars
          matrix.color = "#BB0021FF",# Matrix intersection color
          # text.scale = c(2, 2, 1, 1, 1.2, 3),# text sizes
          # keep.order = TRUE,# preserve input order of sets
          # queries = list(list(query = intersects,
          #                     params = list("Thriller","Action"),
          #                     active = T,color="#EF4143"))
    )
    # Manually export the UpSet figure (A4)
  # Bar chart used in the UpSet figure
    colnames(res)<-c('type','num')
    res$flag<-101:(100+nrow(res))
    gap.barplot(res,y.cols = 'num',btm = 100000,top = 170000,ratio = 2,gap.width = 0.5 , brk.size = 0.4,col='#631879')
    # Manually export to ./graph/upset_bar1.pdf (A4)
  # Another bar chart for dataset counts
    dataset<-colSums(data)
    dataset<-as.data.frame(dataset)
    dataset$flag<-paste0('t',1:nrow(dataset))
    colnames(dataset)<-c('num','flag')
    p<-ggplot(dataset, aes(flag,num)) +
      geom_col(position = position_dodge(width = 0.8),fill='#008280FF',color="black") +
      labs(x = NULL, y = NULL) +
      scale_x_discrete(breaks=dataset$flag,labels=rownames(dataset))+
      theme_bw() + mythem
    # Save
    ggsave(
      filename = './graph/upset_bar2.pdf',
      plot = p,width = 8,height = 6,
      units = 'in'
    )
    remove(data)
  }
}
# Overlap between our identified MNVs and gnomAD
{
  # mnv<-fread('02adjustOrigin/flag_datasets', sep = "\t", quote = '', header = T, stringsAsFactors = F, na.strings = "", data.table = F)
  mymnv<-rowSums(mnv[9:13])
  gnomad<-rowSums(mnv[14:15])
  mylist<-list(mnv=mnv[mymnv!=0,'mnvid'],gnomad=mnv[gnomad!=0,'mnvid'])
  p<-graph_venn(mylist,my_alpha = 0.8)
}
# Density correlation plot
{
  # Normalize to 0~1
  scaleDf<-function(data){
    # Assume a numeric vector input
    # data <- c(1, 2, 3, 4, 5)
    
    # Min/Max
    min_val <- min(data)
    max_val <- max(data)
    
    # Scale to [0,1]
    normalized_data <- (data - min_val) / (max_val - min_val)
    
    # Return scaled values
    return(normalized_data)
  }
  # Load density files
  getdf<-function(url,f){
    mnv<-read.table(paste0(url,'mnv_',f,'_Density.txt'),header = F,sep = '\t',stringsAsFactors = F)
    snv<-read.table(paste0(url,'snv_',f,'_Density.txt'),header = F,sep = '\t',stringsAsFactors = F)
    df<-cbind(snv,mnv$V3)
    colnames(df)<-c('chr','pos','snv','mnv')
    df$flag<-1:nrow(df)
    df$dataset<-f
    df<-df[6:3]
    df$scale_mnv<-scaleDf(df$mnv)
    df$scale_snv<-scaleDf(df$snv)
    return(df)
  }
  # Load all datasets (MHC hg38L:chr6 28510120-33480577)
  {
    df_1000G<-getdf('./03density/03density/','1000G')
    df_GTEx<-getdf('./03density/03density/','GTEx')
  # TCGA MNV missing 78 records on chr18; added manually
    df_TCGA<-getdf('./03density/03density/','TCGA')
    df_UKB20w<-getdf('./03density/03density/','UKB20w')
    df_UKB50w<-getdf('./03density/03density/','UKB50w')
    
    df<-rbind(df_1000G,df_GTEx,df_UKB50w,df_TCGA,df_UKB20w)
    # df$snv<-scaleDf(df$snv)
    # df$mnv<-scaleDf(df$mnv)
  }
  # Correlation
  {
    dfcor<-cor.test(df$snv,df$mnv,method = 'spearm')
    # print(c(dfcor$estimate,dfcor$p.value))
    df$dataset<-factor(df$dataset)
    p<-ggplot(data=df,aes(x=snv,y=mnv,color=dataset))+
      geom_point(pch=16,size=3,alpha=0.5)+
      labs(title = paste0('R:',round(dfcor$estimate, 2)),x='Density of SNVs (SNVs/1MB)',y='Density of MNVs (MNVs/1MB)')+ 
      scale_color_manual(values = c("#7E277B", "#FCA41C","#0D8040","#3353A3","#EF251F"))+
      theme_bw() + mythem
    # Save
    ggsave(
      filename = './graph/cor.pdf',
      plot = p,width = 6,height = 5,
      units = 'in'
    )
    # Optionally, draw a legend-free version for figure assembly
  }
}
# Why only identify up to 5-joint MNVs
{
  ## Plot with error bars (read length = 100bp)
  data<-read.table('04why5joint/statistic.file5',header = F,stringsAsFactors = F,sep = '\t')
  colnames(data)<-c('group','sample','value')
  data<-data %>%
    group_by(group) %>%
    summarise(mean = mean(value), 
              se = sqrt(var(value)/length(value)) ,
              .groups = "drop")  %>% 
    mutate(group = as.character(group)) 
  p<-ggplot(data,aes(factor(group),mean)) +
    scale_y_continuous(breaks=seq(0.4, 1, by = 0.1),limits = c(0.4, 1)) +
    geom_point(size=3,aes(color=group))+
    geom_errorbar(aes(ymin = mean-se,ymax =mean+se,
                      group = group,color=group),
                  width = 0.1)+
    geom_hline(yintercept = c(0.5), linetype = "dashed", color = c("red"))+
    scale_color_nejm()+
    labs(x='MNV Type (N joint MNV) (The Read Length = 100bp)',y='Mean Mapping Rate (%)')+theme_bw() +mythem
  ggsave(
    filename = 'graph/select5point100bp.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  ## Plot with error bars (read length = 150bp)
  data<-read.table('04why5joint/statistic.file',header = F,stringsAsFactors = F,sep = '\t')
  colnames(data)<-c('group','sample','value')
  data<-data %>%
    group_by(group) %>%
    summarise(mean = mean(value), 
              se = sqrt(var(value)/length(value)) ,
              .groups = "drop")  %>% 
    mutate(group = as.character(group))
  p<-ggplot(data,aes(factor(group),mean)) +
    scale_y_continuous(breaks=seq(0.5, 1, by = 0.1),limits = c(0.5, 1)) +
    geom_point(size=3,aes(color=group))+
    geom_errorbar(aes(ymin = mean-se,ymax =mean+se,
                      group = group,color=group),
                  width = 0.1)+
    scale_color_nejm()+
    labs(x='MNV Type (N joint MNV) (The Read Length = 150bp)',y='Mean Mapping Rate (%)')+theme_bw() +mythem
  ggsave(
    filename = 'graph/select5point150bp.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
# AC error correction: how many MNVs corrected and how many low-joint miscalls removed
{
  data<-read.table('05ACError/merge.txt',header = F,stringsAsFactors = F,sep = '\t',row.names = 1)
  data<-rbind(data,colSums(data))
  data$f1<-data$V3/data$V2
  data$f2<-data$V4/data$V2
  data$f3<-data$V4/data$V3
  data<-rbind(data,c(0,0,0,median(data[1:37,'f1']),median(data[1:37,'f2']),median(data[1:37,'f3'])))
  colnames(data)<-c('total','adjust','filter','a/t','f/t','f/a')
  write.table(data,'graph/ACError.txt',sep = '\t',quote = F)
  # On average, corrected about 14% of MNV frequency
}
