library(ggplot2)
library(ggsci)
library(scales)
library(data.table)
library(tidyverse)
library(venn)
library(UpSetR)
# Extract AAAS color palette
color1 = pal_aaas()(10)
show_col(color1)
colorPalette<-c("#e41a1c","#377eb8","#4daf4a","9ecae1","#6baed6","#4292c6")
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
}
# Load packages and utilities for drawing truncated (gap) barplots
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
          ## Auto-calculate btm/top if not provided
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
  ## Function to draw error bars
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
  ## Function to automatically compute break positions
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
# Read data
{
  readDf<-function(x){
    df<-read.table(paste0(x,'.mnv'),header = F,stringsAsFactors = F,sep = '\t')
    df<-df[1:8]
    colnames(df)<-c('chr','pos','mnvid','ref','alt','rsid','distance','mnvtype')
    df<-df[!duplicated(df),]
    df$flag<-apply(df[c(2,4,5)],1,paste,collapse='-')
    df$origin<-x
    return(df)
  }
  bcftools<-readDf('bcftools')
  MAC<-readDf('MAC')
  MACARON<-readDf('MACARON')
  MNVAnno<-readDf('MNVAnno')
  Qing<-readDf('Qing')
  vardict<-readDf('vardict')
}
# Count summary
{
  # Create example vector
  mnv_flag <- c("bcftools", "MACARON", "MAC","MNVAnno","Qing","vardict")
  mnv_num <- c(nrow(bcftools), nrow(MACARON), nrow(MAC),nrow(MNVAnno),nrow(Qing),nrow(vardict))
  # Plot
  p<-ggplot(my_df, aes(x = factor(flag, levels = mnv_flag), y = num)) +
    coord_flip()+
    # geom_col(position = position_dodge(width = 0.8), fill = '#008280FF', color = "black") +
    geom_bar(stat = "identity", fill = "#008280FF", width = 0.6, color = "black") +
    scale_y_continuous(breaks=seq(0, 16000, by = 4000),limits = c(0, 16000)) +
    geom_text(aes(label = num), vjust = -0.5,angle=-90, color = "black")+
    labs(x = NULL, y = NULL) +
    scale_x_discrete(limits = mnv_flag) +
    theme_bw() + mythem
  # gap.barplot(my_df,y.cols = 'num',btm = 20,top = 1000,ratio = 2,gap.width = 0.5 , brk.size = 0.4,col='#631879')
  ggsave(
    filename = '../graph/01overlap_1.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
# overlap
{
  # venn
  {
    # Exclude bcftools and MACARON from the overlap due to too few calls
    x <- list(
              MAC=MAC$flag,
              MNVAnno=MNVAnno$flag,
              Qing=Qing$flag,
              vardict=vardict$flag
    )
    p<-venn::venn(x,zcolor=2:5,box=T)
  }
  ## upset
  {
    p2<-p$counts[2:16]
    names(p2)<-gsub(':','&',rownames(p)[2:16])
    data2 <- fromExpression(p2)
    upset(data2, nsets = 7,  sets =  c('MAC','MNVAnno','Qing','vardict'),
          matrix.color ="#b35806", 
          # main.bar.color = colorPalette[1:4],
          # sets.bar.color = c("#e41a1c","#377eb8","#4daf4a"),
          order.by = c("freq"), decreasing = c(TRUE),
          point.size = 4,  line.size = 1.3,  mainbar.y.label = "IntersectionSize", 
          keep.order = TRUE,
          sets.x.label = "", mb.ratio = c(0.60, 0.40), text.scale = c(2, 2, 0.5, 0.5,2, 3))
    # Plot
    data2<-data.frame(flag = gsub(':','&',rownames(p)[2:16]), num = p$counts[2:16])
    data2 <- data2[order(-data2$num), ]
    data2<-data2[data2$num>0,]
    p<-ggplot(data2, aes(x = factor(flag, levels = flag), y = num)) +
      geom_bar(stat = "identity", fill = "#5F2D73", width = 0.6, color = "black") +
      scale_y_continuous(breaks=seq(0, 15000, by = 3000),limits = c(0, 15000)) +
      geom_text(aes(label = num), vjust = -0.5, color = "black")+
      labs(x = NULL, y = NULL) +
      # scale_x_discrete(limits = mnv_flag) +
      theme_bw() + mythem
    ggsave(
      filename = '../graph/01overlap_2.pdf',
      plot = p,width = 8,height = 6,
      units = 'in'
    )
  }
}
# Time
{
  library(tidyverse)
  library(ggsci)
  df_time<-read.table('time.txt',header = T,stringsAsFactors = F,sep = '\t')
  df_time[1:15,'time']<- df_time[1:15,'time']+120 # To plot on the same figure
  ## Plot with error bars
  df_time<-df_time %>%
      group_by(method) %>%
      summarise(mean = mean(time), 
                se = sqrt(var(time)/length(time)) ,
                .groups = "drop")
  p<-ggplot(df_time,aes(factor(method, levels = c('BCFtools','MNVAnno','Qing','MACARON','MAC','varDic')),mean)) +
    scale_y_continuous(breaks=seq(0, 250, by = 50),limits = c(0, 250)) +
    # geom_line(aes(linetype=method,group = method,color=method))+
    geom_point(size=3,aes(color=method))+
    # geom_point(shape=21, fill="#69b3a2", size=3,aes(color=tools))+
    geom_errorbar(aes(ymin = mean-se,ymax =mean+se,
                      group = method,color=method),
                  width = 0.1)+
    scale_color_nejm()+
    labs(x='',y='Time Costs')+theme_bw() +mythem
  ggsave(
    filename = '03time_bar.pdf',
    plot = p,width = 12,height = 5,
    units = 'in'
  )
}
# Memory
{
  df_mem<-read.table('mem.txt',header = T,stringsAsFactors = F,sep = '\t')
  df_mem[1:10,'mem']<- df_mem[1:10,'mem']*100+2000 # To plot on the same figure
  ## Plot with error bars
  df_mem<-df_mem %>%
    group_by(method) %>%
    summarise(mean = mean(mem), 
              se = sqrt(var(mem)/length(mem)) ,
              .groups = "drop") 
  p<-ggplot(df_mem,aes(factor(method, levels = c('MNVAnno','varDic','BCFtools','Qing','MACARON','MAC')),mean)) +
    scale_y_continuous(breaks=seq(0, 3600, by = 600),limits = c(0, 3600)) +
    # geom_line(aes(linetype=method,group = method,color=method))+
    geom_point(size=3,aes(color=method))+
    # geom_point(shape=21, fill="#69b3a2", size=3,aes(color=tools))+
    geom_errorbar(aes(ymin = mean-se,ymax =mean+se,
                      group = method,color=method),
                  width = 0.1)+
    geom_hline(yintercept = c(26, 714), linetype = "dashed", color = c("blue","red"))+
    scale_color_nejm()+
    labs(x='',y='Memory Costs')+theme_bw() +mythem
  ggsave(
    filename = '04mem_bar.pdf',
    plot = p,width = 12,height = 5,
    units = 'in'
  )
}
# Accuracy, sensitivity, specificity, ROC (use 10 reads; only the four VCF-based methods)
{
  # Read data
  readDf<-function(x){
    df<-read.table(paste0(x,'_SN.txt'),header = T,stringsAsFactors = F,sep = '\t')
    colnames(df)[1]<-x
    return(df)
  }
  # Compute accuracy, precision, sensitivity, specificity
  getThreeVal<-function(data,aa_results){
    # Create dataset (example)
    # data <- data.frame(
    #   id = c("mnv1", "mnv2", "mnv3", "mnv4", "mnv5", "mnv6", "mnv7", "mnv8", "mnv9", "mnv10", "mnv11"),
    #   isture = c(1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0)
    # )
    # 
    # Calls identified by method aa (example)
    # aa_results <- c("mnv2", "mnv3", "mnv6", "mnv7", "mnv8")
    
    # Calculate TP (True Positive)
    TP <- sum(data$id %in% aa_results & data$isture == 1)
    
    # Calculate TN (True Negative)
    TN <- sum(!(data$id %in% aa_results) & data$isture == 0)
    
    # Calculate FP (False Positive)
    FP <- sum(data$id %in% aa_results & data$isture == 0)
    
    # Calculate FN (False Negative)
    FN <- sum(!(data$id %in% aa_results) & data$isture == 1)
    
    # Calculate precision
    precision <- TP / (TP + FP)
    
    # Calculate accuracy
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    
    # Calculate sensitivity
    sensitivity <- TP / (TP + FN)
    
    # Calculate specificity
    specificity <- TN / (TN + FP)
    
    # Output (examples)
    # cat("Accuracy:", accuracy, "\n")
    # cat("Sensitivity:", sensitivity, "\n")
    # cat("Specificity:", specificity, "\n")
    return(c(accuracy,precision,sensitivity,specificity))
  }
  # Merge results from the four methods
  getThreeMerge<-function(df2,cutoff){
    data<-df2[c(1,7)]
    colnames(data)<-c('id','isture')
    aa_results<-df2[df2$MAC>=cutoff,1]
    effect_MAC<-getThreeVal(data,aa_results)
    aa_results<-df2[df2$MNVAnno>=cutoff,1]
    effect_MNVAnno<-getThreeVal(data,aa_results)
    aa_results<-df2[df2$Qing>=cutoff,1]
    effect_Qing<-getThreeVal(data,aa_results)
    aa_results<-df2[df2$vardict>=cutoff,1]
    effect_vardict<-getThreeVal(data,aa_results)
    # Use data.frame() to create a data frame
    my_flag <- c("MAC","MNVAnno","Qing","vardict")
    my_cutoff <- c(rep(cutoff,4))
    my_accuracy <- c(effect_MAC[1],effect_MNVAnno[1],effect_Qing[1],effect_vardict[1])
    my_precision<- c(effect_MAC[2],effect_MNVAnno[2],effect_Qing[2],effect_vardict[2])
    my_sensitivity <- c(effect_MAC[3],effect_MNVAnno[3],effect_Qing[3],effect_vardict[3])
    my_specificity <- c(effect_MAC[4],effect_MNVAnno[4],effect_Qing[4],effect_vardict[4])
    my_df <- data.frame(method = my_flag, cutoff = my_cutoff, accuracy = my_accuracy, precision=my_precision , sensitivity = my_sensitivity, specificity = my_specificity)
    return(my_df)
  }
  # Load dataset
  df<-read.table('./check/merge/res_adjust_cutoff.txt',header = T,stringsAsFactors = F,sep = '\t')
  df$flag<-apply(df[2:4],1,paste,collapse='-')
  df<-df[c(1:5,9,11)] # cut10
  colnames(df)[6]<-'cut10'
  df[df$cut10>=1,'cut10']<-1
  df<-merge(df,readDf('MAC'),by='flag',all.x = T,sort=F)
  df<-merge(df,readDf('MNVAnno'),by='flag',all.x = T,sort=F)
  df<-merge(df,readDf('Qing'),by='flag',all.x = T,sort=F)
  df<-merge(df,readDf('vardict'),by='flag',all.x = T,sort=F)
  # Consider only MNV present in VCF
  ref<-read.table('check/total_filter.mnv',header = F,stringsAsFactors = F,sep = '\t')
  ref$flag<-apply(ref[2:4],1,paste,collapse='-')
  res<-merge(ref[c(5,1)],df,by='flag',all.x = T,sort=F)
  df<-res[c(1,3:12)]
  df[is.na(df)]<-0
  # Get metrics
  my_df<-getThreeMerge(df,1)
  for (i in 2:25) {
    my_df<-rbind(my_df,getThreeMerge(df,i))
  }
  my_df$Mspecificity<-1-my_df$specificity
  # Line plots
  ggplot(my_df, aes(x = cutoff , y = accuracy, color = method)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = 1:25)+
    xlab("cutoff") +
    ylab("accuracy") +
    ggtitle("Group Line Chart")+theme_bw()+mythem
  ggplot(my_df, aes(x = cutoff , y = sensitivity, color = method)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = 1:25)+
    xlab("cutoff") +
    ylab("sensitivity") +
    ggtitle("Group Line Chart")+theme_bw()+mythem
  ggplot(my_df, aes(x = cutoff , y = specificity, color = method)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = 1:25)+
    xlab("cutoff") +
    ylab("specificity") +
    ggtitle("Group Line Chart")+theme_bw()+mythem
  # Use cutoff = 1 for plotting (i.e., at least one supporting sample)
  my_df<-my_df[1:4,]
  p<-ggplot(my_df, aes(x=factor(method, levels = c('MNVAnno','Qing','vardict','MAC')), y=accuracy)) +
    geom_point(size=3,aes(color=method))+
    scale_y_continuous(breaks=seq(0, 1, by = 0.2),limits = c(0, 1)) +
    labs(x = NULL, y = 'accuracy') +
    theme_bw()+mythem
  ggsave(
    filename = '../graph/05accuracy.pdf',
    plot = p,width = 12,height = 5,
    units = 'in'
  )
  p<-ggplot(my_df, aes(x=factor(method, levels = c('MNVAnno','Qing','vardict','MAC')), y=sensitivity)) +
    geom_point(size=3,aes(color=method))+
    scale_y_continuous(breaks=seq(0, 1, by = 0.2),limits = c(0, 1)) +
    labs(x = NULL, y = 'sensitivity') +
    theme_bw()+mythem
  ggsave(
    filename = '../graph/05sensitivity.pdf',
    plot = p,width = 12,height = 5,
    units = 'in'
  )
  p<-ggplot(my_df, aes(x=factor(method, levels = c('vardict','MAC','MNVAnno','Qing')), y=specificity)) +
    geom_point(size=3,aes(color=method))+
    scale_y_continuous(breaks=seq(0, 1, by = 0.2),limits = c(0, 1)) +
    labs(x = NULL, y = 'specificity') +
    theme_bw()+mythem
  ggsave(
    filename = '../graph/05specificity.pdf',
    plot = p,width = 12,height = 5,
    units = 'in'
  )
}

