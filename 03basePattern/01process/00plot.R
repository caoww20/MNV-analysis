library(ggplot2)
library(ggsci)
library(scales)

# Extract AAAS palette & themes
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
  mythem_angel<-theme(panel.grid=element_blank(),
                plot.margin=unit(rep(2,4),'lines'),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                title = element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.text.x=element_text(vjust=0,size=16,face = "bold", color='black',angle = 45),
                axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black')
  )
  mythem2<-theme(panel.grid=element_blank(),
                plot.margin=unit(rep(2,4),'lines'),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                title = element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.text.x=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.x=element_text(vjust=0, size=16,face = "bold", color='black'),
                axis.text.y=element_text(vjust=1,size=16,face = "bold", color='black'),
                axis.title.y=element_text(vjust=2, size=16,face = "bold", color='black')
  )
}
# motif utilities
{
  # install.packages('ggseqlogo')
  library('ggseqlogo')
  
  letterMatrix <- function(input){
    # Ensure kmers are the same length characters (ggseqlogo)
    # All sequences must have identical length
    seq.len = sapply(input, nchar) # length per sequence
    num_pos = seq.len[1]
    if(! all(seq.len == num_pos)) {
      stop('Sequences in alignment must have identical lengths')
    }
    
    # Construct matrix of letters (ggseqlogo)
    split = unlist( sapply(input, function(seq){strsplit(seq, '')}) )
    
    t( matrix(split, seq.len, length(split)/num_pos) )
  }
  make_ppm <- function(seqs, ppm=TRUE, seq_type="dna") {
    # seqs: vector of DNA/RNA sequences
    # ppm: return PPM (default) else PFM
    # seq_type: "dna" or "rna"
    
    letter_mat = letterMatrix(seqs)
    
    # Get namespace (ggseqlogo)
    if(seq_type == "dna") {
      namespace = c("A", "T", "G", "C") 
    } else if (seq_type == "rna" ) {
      namespace = c("A", "U", "G", "C") 
    } else {
      stop('Wrong seq_type! Must be one of "dna" and "rna".')
    }
    
    # Construct PWM (ggseqlogo)
    pfm_mat = apply(letter_mat, 2, function(pos.data){
      # Get base counts per position
      t = table(pos.data)
      # Match to namespace order
      ind = match(namespace, names(t))
      # Reordered column
      col = t[ind]
      col[is.na(col)] = 0
      names(col) = namespace
      
      if(ppm) { # Convert to PPM
        col = col / sum(col)      
      }
      col
    })
    
    num_pos = nchar(seqs[1])
    colnames(pfm_mat) = 1:num_pos
    pfm_mat
    
  }
  pfm2ppm <- function(pfm) {
    ppm <- apply(pfm, 2, function(col) {col / sum(col)} )
    return(ppm)
  }
  ppm2pwm <- function(ppm) {
    pwm <- log2(ppm / 0.25)
    pwm[is.infinite(pwm)] <- 0 # Replace -Inf with 0 when freq=0
    return(pwm)
  }

}
# Base complement/transformation helpers
{
  changeAT<-function(x){
    if (x=='A') {
      x='T'
    }else if (x=='T'){
      x='A'
    }else if (x=='C'){
      x='G'
    }else if (x=='G'){
      x='C'
    }
    return(x)
  }
  transATCG<-function(x){
    # x='A,T'
    x=str_split(x,',',simplify = T)
    m=c()
    for (i in length(x):1) {
      m=c(m,changeAT(x[i]))
    }
    x=paste(m,collapse = ',')
    return(x)
  }
  transATCG2<-function(x){
    # x='A,C->T,T'
    x=str_split(x,'->',simplify = T)
    a=str_split(x[1],',',simplify = T)
    m=c()
    for (i in length(a):1) {
      m=c(m,changeAT(a[i]))
    }
    a=paste(m,collapse = ',')
    
    b=str_split(x[2],',',simplify = T)
    m=c()
    for (i in length(b):1) {
      m=c(m,changeAT(b[i]))
    }
    b=paste(m,collapse = ',')
    
    x=paste0(a,'->',b)
    return(x)
  }
  getPattern<-function(x){
  a=c(x[1],x[2]) # a=A,A->C,C; b=T,T->G,G
    a=sort(a)
    return(a[1])
  }
  getTrans<-function(x){
  # AAAC -> GTTT
    x=unlist(strsplit(x, ""))
    m=c()
    for (i in length(x):1) {
      m=c(m,changeAT(x[i]))
    }
    x=paste(m,collapse = '')
    return(x)
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
# Load base data
{
  data<-fread('mnv_5_adjust.txt', sep = "\t", quote = '', header = T, stringsAsFactors = F, na.strings = "", data.table = F)
}
# Global summary
{
  res<-table(data$distance,data$mnvtype)
  res<-as.data.frame(res,stringsAsFactors = F)
  colnames(res)<-c('distance','mnvtype','number') 
  res$distance<-as.numeric(res$distance)
  res$mnvtype<-as.numeric(res$mnvtype)
  write.table(res,'../graph/summary1.txt',quote = F,row.names = F,sep = '\t')
  # Convert summary1.txt to 1.tsv
  tmp<-cbind(res[res$mnvtype==2,c(1,3)],res[res$mnvtype==3,3],res[res$mnvtype==4,3],res[res$mnvtype==5,3])
  colnames(tmp)<-c('distance','2','3','4','5')
  write.table(tmp,'../graph/1.tsv',quote = F,row.names = F,sep = '\t')
  # Manual tweaks to avoid sharp spikes
  res[res$distance==1 & res$mnvtype==3,'number']<-5000
  res[res$distance==2 & res$mnvtype==4,'number']<-300
  res[res$distance==3 & res$mnvtype==5,'number']<-100
  # Value scaling
  res$number<-log10(res$number+1)
  
  # Curves per MNV type
  line2=res[res$mnvtype==2,]
  line3=res[res$mnvtype==3,]
  line4=res[res$mnvtype==4,]
  line5=res[res$mnvtype==5,]
  point2 = spline(line2$distance,line2$number,1000)
  point3 = spline(line3$distance,line3$number,1000)
  point4 = spline(line4$distance,line4$number,1000)
  point5 = spline(line5$distance,line5$number,1000)
  
  abab=data.frame(point2$x,point2$y,point3$y,point4$y,point5$y)
  colnames(abab)<-c('distance',paste0('j',2:5))
  abab[abab$distance<2,'j3']<-0
  abab[abab$distance<3,'j4']<-0
  abab[abab$distance<4,'j5']<-0
  
  
  p<-ggplot(abab,aes(x = distance)) +
    geom_point(aes(y = j2),col='#B0282E')+
    geom_point(aes(y = j3),col='#FCA41C')+
    geom_point(aes(y = j4),col='#7E277B')+
    geom_point(aes(y = j5),col='#3353A3')+
    scale_x_continuous(breaks=c(1:10))+
    ylim(c(2,6))+theme_bw()+mythem
  p
  ggsave(
    filename = '../graph/figA.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
## Correlation with/without adjacent constraints
{
  # Draw correlation
  corrdf<-read.table('../graph/1.tsv',sep = '\t',header = T,stringsAsFactors = T)
  
  X2<-corrdf$X2 # max=2 (excluding near-end points)
  X3<-corrdf$X3 # max=4
  X4<-corrdf$X4 # max=6
  X5<-corrdf$X5 # max=8
  
  # Revised version
  getCor<-function(X2,X3,X4,X5,flag){
    a=cor(X2[1:8],X3[3:10],method = flag)
    a1=cor(X2[1:6],X4[5:10],method = flag)
    a2=cor(X2[1:4],X5[7:10],method = flag)
    
    b1=cor(X3[2:8],X4[4:10],method = flag)
    b2=cor(X3[2:6],X5[6:10],method = flag)
    
    c2=cor(X4[3:8],X5[5:10],method = flag)
    
    res<-matrix(c(1,a,a1,a2,a,1,b1,b2,a1,b1,1,c2,a2,b2,c2,1),nrow = 4,ncol = 4,byrow = F)
    colnames(res)<-c('X2','X3','X4','X5')
    rownames(res)<-c('X2','X3','X4','X5')
    return(res)
  }
  
  res<-getCor(X2,X3,X4,X5,'pearson')
  # res<-getCor(X2,X3,X4,X5,'spearman')
  library(corrplot)
  corrplot(res)
  # Manually export to corr.pdf (size 640x560)
}
# Adjacent 2-joint ref->alt matrix
{
  joint2<-data[data$mnvtype==2 & data$distance==1,]
  changeMatrix<-table(joint2[2:3])
  changelist<-as.data.frame(changeMatrix)
  changelist<-changelist[changelist$Freq>0,]
  changelist$trans_refs<-apply(changelist['refs'],1,transATCG)
  changelist$trans_alts<-apply(changelist['alts'],1,transATCG)
  changelist$type<-paste0(changelist$refs,'->',changelist$alts)
  changelist$type_trans<-paste0(changelist$trans_refs,'->',changelist$trans_alts)
  changelist<-changelist[c(6,7,3)]
  reslist<-list()
  for (i in 1:nrow(changelist)) {
    if (changelist[i,1] != 'F') {
      type<-changelist[i,1]
      type_trans<-changelist[i,2]
      a<-changelist[changelist$type==type,3]
      b<-changelist[changelist$type_trans==type_trans,3]
      reslist[[type]]<-c(type,a+b)
      changelist[changelist$type==type_trans,1]<-'F'
    }
  }
  changelist<-as.data.frame(do.call(rbind,reslist))
  colnames(changelist)<-c('type','number')
  changelist$number<-as.numeric(changelist$number)
  changelist<-changelist[order(changelist$number),]
  write.table(changelist,'../graph/summary2.txt',sep = '\t',row.names = F,quote = F)

  ## plot
  patterndf<-rbind(changelist[1:7,],changelist[(nrow(changelist)-5):nrow(changelist),])
  patterndf$log10Value<-log10(patterndf$number+1)
  patterndf$id<-1:nrow(patterndf)
  patterndf[7,1]<-'...'
  patterndf[7,2]<-0
  
  gap.barplot(patterndf,y.cols = 'number',btm = 1500,top = 15000,ratio = 4,gap.width = 0.5 , brk.size = 0.4)
  # Manually export to figB.pdf (8x6)
}
# Non-adjacent 2-joint ref->alt matrix
{
  joint2<-data[data$mnvtype==2 & data$distance!=1,]
  changeMatrix<-table(joint2[2:3])
  changelist<-as.data.frame(changeMatrix)
  changelist<-changelist[changelist$Freq>0,]
  changelist$trans_refs<-apply(changelist['refs'],1,transATCG)
  changelist$trans_alts<-apply(changelist['alts'],1,transATCG)
  changelist$type<-paste0(changelist$refs,'->',changelist$alts)
  changelist$type_trans<-paste0(changelist$trans_refs,'->',changelist$trans_alts)
  changelist<-changelist[c(6,7,3)]
  reslist<-list()
  for (i in 1:nrow(changelist)) {
    if (changelist[i,1] != 'F') {
      type<-changelist[i,1]
      type_trans<-changelist[i,2]
      a<-changelist[changelist$type==type,3]
      b<-changelist[changelist$type_trans==type_trans,3]
      reslist[[type]]<-c(type,a+b)
      changelist[changelist$type==type_trans,1]<-'F'
    }
  }
  changelist<-as.data.frame(do.call(rbind,reslist))
  colnames(changelist)<-c('type','number')
  changelist$number<-as.numeric(changelist$number)
  changelist<-changelist[order(changelist$number),]
  ## plot
  patterndf<-rbind(changelist[1:7,],changelist[(nrow(changelist)-5):nrow(changelist),])
  patterndf$id<-1:nrow(patterndf)
  patterndf[7,1]<-'...'
  patterndf[7,2]<-0
  
  gap.barplot(patterndf,y.cols = 'number',btm = 10000,top = 100000,ratio = 4,gap.width = 0.5 , brk.size = 0.4)
  # Manually export to figB_sub.pdf (8x6)
}
# Proportion of mechanisms per joint type
{
  mecPer<-list()
  for (i in 2:5) {
    mecPer[[i-1]]<-colSums(data[data$mnvtype==i,12:14])
  }
  mecPer<-do.call(rbind,mecPer)
  mecPer<-as.data.frame(mecPer)
  mecPer<-mecPer/rowSums(mecPer)
  mecPer$mnvtype<-paste0(2:5)
  
  mecPer<-reshape2::melt(mecPer,id=c('mnvtype'))
  colnames(mecPer)<-c('mnvtype','mechanism','value')
  p<-ggplot(mecPer,aes(x=mnvtype,y=100*value, fill = mechanism))+
  geom_col(position="stack",width = 0.5)+ # stacked bar
    labs(x='MNV Type (N joint MNV)',y='Relative Abundance(%)')+ 
    scale_fill_manual(values = c("snv_event" = "#2A5888", "one_step" = "#4396B1", "repeat" = "#89CEED"))+
    theme_bw()+mythem
  ggsave(
    filename = '../graph/figC.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
# SNV event analysis
{
  # Note: length(flag)>0 is fine for human but may break for animals with many zeros
  getPatternFreq<-function(df,cutoff){
    # res=list()
    df$pattern<-paste0(df$refs,'->',df$alts)
    df<-df[c(15,12,10,11)]
    # No collapsing yet
    df_matrix<-as.data.frame(table(df[1:2]),stringsAsFactors = F)
    df_matrix<-cbind(df_matrix[df_matrix$snv_event==0,c(1,3)],df_matrix[df_matrix$snv_event==1,3])
    colnames(df_matrix)<-c('pattern','notrepeat','repeat')
    df_matrix<-rbind(df_matrix,c(0,sum(df_matrix$notrepeat),sum(df_matrix$`repeat`)))
    df_matrix[nrow(df_matrix),1]<-'all'
    df_matrix$freq<-df_matrix$`repeat`/(df_matrix$`repeat`+df_matrix$notrepeat)
    # res[[1]]=df_matrix[order(-df_matrix$freq),] # optional
    # Collapse reverse-complement patterns
    df_matrix<-df_matrix[1:(nrow(df_matrix)-1),]
    df_matrix$pattern2<-apply(df_matrix[1], 1, transATCG2)
    df_matrix<-df_matrix[c(1,5,2:4)]
    df_matrix$onePattern<-apply(df_matrix[1:2],1,getPattern)
    df_matrix<-df_matrix[c(6,3,4)]
    df_list=list()
    for (i in unique(df_matrix$onePattern)) {
      df_list[[i]]<-colSums(df_matrix[df_matrix$onePattern==i,2:3])
    }
    df_matrix<-as.data.frame(do.call(rbind,df_list))
    df_matrix$pattern<-rownames(df_matrix)
    df_matrix<-df_matrix[c(3,1,2)]
    df_matrix<-rbind(df_matrix,c(0,sum(df_matrix$notrepeat),sum(df_matrix$`repeat`)))
    df_matrix[nrow(df_matrix),1]<-'all'
    df_matrix$total<-df_matrix$`repeat`+df_matrix$notrepeat
    df_matrix$freq<-df_matrix$`repeat`/df_matrix$total
    df_matrix<-df_matrix[order(-df_matrix$total),]
    flag=which(df_matrix$`repeat`==0)
    if (length(flag)>0) {
      df_matrix<-df_matrix[1:(which(df_matrix$`repeat`==0)[1]-1),]
    }
    df_matrix$realFreq<-df_matrix$freq/df_matrix[df_matrix$pattern=='all','freq']
    df_matrix<-df_matrix[df_matrix$total>cutoff,]
    df_matrix<-df_matrix[order(-df_matrix$freq),]
    return(df_matrix)
  }
  # Without restricting to contiguous mutations; run twice if needed (total vs MAF 0.05)
  m2<-data[data$mnvtype==2,] 
  m2_res<-getPatternFreq(m2,100)
  m3<-data[data$mnvtype==3,] 
  m3_res<-getPatternFreq(m3,100)
  m4<-data[data$mnvtype==4 ,] 
  m4_res<-getPatternFreq(m4,100)
  m5<-data[data$mnvtype==5,] 
  m5_res<-getPatternFreq(m5,100)
  
  # Scatter plots
  getTop5<-function(df){
    point<-rev(df[rev(c(1:5,(nrow(df)-5):nrow(df))),])
    point$id<-1:nrow(point)
    return(point)
  }
  pFreq<-function(df){
    point<-getTop5(df)
    backgroudValue<-df[df$pattern=='all','freq']
    p<-ggplot(point,aes(x = id,y = freq)) +
      geom_point(aes(size=total,color=freq))+
      # scale_colour_gradientn(colors = c("#FFD133", "#FF8D33", "#B0282E"), values=c(0,0.5,1), limits=c(0,1), guide="colourbar")
      scale_colour_gradientn(colors = c("#FFD133", "#FF8D33", "#B0282E"))+ # color gradient
      scale_size_continuous(range=c(3,5))+ # point size range
      scale_x_continuous(breaks=c(1:nrow(point)),labels=point$pattern)+
      # scale_y_continuous(breaks=c(0:1))
      ylim(c(0,1))+
      geom_hline(aes(yintercept=backgroudValue),linetype=5,col="red")+
      # labs(x='MNV Pattern',y='Fraction of MNVs in repetitive contexts')
      labs(x='',y='')+
      theme_bw()+mythem2
    return(p)
  }
  p<-pFreq(m2_res)
  ggsave(
    filename = '../graph/snv_event/2joint_freq.pdf',
    plot = p,width = 8,height = 3,
    units = 'in'
  )
}
# SNV event (adjacent MNV only)
{
  m2<-data[data$mnvtype==2 & data$distance==1,] 
  m2_res<-getPatternFreq(m2,100)
  # Scatter plot
  p<-pFreq(m2_res)
  ggsave(
    filename = '../graph/snv_event/adjacent_2joint_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
# One-step mechanism
{
  getPatternFreq<-function(df,cutoff){
    # res=list()
    df$pattern<-paste0(df$refs,'->',df$alts)
    df<-df[c(15,13,10,11)]
    # No collapsing yet
    df_matrix<-as.data.frame(table(df[1:2]),stringsAsFactors = F)
    df_matrix<-cbind(df_matrix[df_matrix$one_step==0,c(1,3)],df_matrix[df_matrix$one_step==1,3])
    colnames(df_matrix)<-c('pattern','notrepeat','repeat')
    df_matrix<-rbind(df_matrix,c(0,sum(df_matrix$notrepeat),sum(df_matrix$`repeat`)))
    df_matrix[nrow(df_matrix),1]<-'all'
    df_matrix$freq<-df_matrix$`repeat`/(df_matrix$`repeat`+df_matrix$notrepeat)
    # res[[1]]=df_matrix[order(-df_matrix$freq),] # optional
    # Collapse reverse-complement patterns
    df_matrix<-df_matrix[1:(nrow(df_matrix)-1),]
    df_matrix$pattern2<-apply(df_matrix[1], 1, transATCG2)
    df_matrix<-df_matrix[c(1,5,2:4)]
    df_matrix$onePattern<-apply(df_matrix[1:2],1,getPattern)
    df_matrix<-df_matrix[c(6,3,4)]
    df_list=list()
    for (i in unique(df_matrix$onePattern)) {
      df_list[[i]]<-colSums(df_matrix[df_matrix$onePattern==i,2:3])
    }
    df_matrix<-as.data.frame(do.call(rbind,df_list))
    df_matrix$pattern<-rownames(df_matrix)
    df_matrix<-df_matrix[c(3,1,2)]
    df_matrix<-rbind(df_matrix,c(0,sum(df_matrix$notrepeat),sum(df_matrix$`repeat`)))
    df_matrix[nrow(df_matrix),1]<-'all'
    df_matrix$total<-df_matrix$`repeat`+df_matrix$notrepeat
    df_matrix$freq<-df_matrix$`repeat`/df_matrix$total
    df_matrix<-df_matrix[order(-df_matrix$total),]
    flag=which(df_matrix$`repeat`==0)
    if (length(flag)>0) {
      df_matrix<-df_matrix[1:(which(df_matrix$`repeat`==0)[1]-1),]
    }
    df_matrix$realFreq<-df_matrix$freq/df_matrix[df_matrix$pattern=='all','freq']
    df_matrix<-df_matrix[df_matrix$total>cutoff,]
    df_matrix<-df_matrix[order(-df_matrix$freq),]
    return(df_matrix)
  }
  # Not restricting to contiguous MNVs: run twice (total vs MAF 0.05)
  m2<-data[data$mnvtype==2,] 
  m2_res<-getPatternFreq(m2,100)
  m3<-data[data$mnvtype==3,] 
  m3_res<-getPatternFreq(m3,100)
  m4<-data[data$mnvtype==4 ,] 
  m4_res<-getPatternFreq(m4,100)
  m5<-data[data$mnvtype==5,] 
  m5_res<-getPatternFreq(m5,100)
  
  # Draw scatter plot
  getTop5<-function(df){
    point<-rev(df[rev(c(1:5,(nrow(df)-5):nrow(df))),])
    point$id<-1:nrow(point)
    return(point)
  }
  pFreq<-function(df){
    point<-getTop5(df)
    backgroudValue<-df[df$pattern=='all','freq']
    p<-ggplot(point,aes(x = id,y = freq)) +
      geom_point(aes(size=total,color=freq))+
      scale_colour_gradientn(colors = c("#FFD133", "#FF8D33", "#B0282E"))+ # 绘制颜色渐变
      scale_size_continuous(range=c(3,5))+ # 将点大小限制在某一个范围内 连续型
      # scale_size_continuous(range=c(3,5))+ # 将点大小限制在某一个范围内 离散型
      scale_x_continuous(breaks=c(1:nrow(point)),labels=point$pattern)+
      # scale_y_continuous(breaks=c(0:1))+
      # ylim(c(ymin,1))+
      ylim(c(0,1))+
      geom_hline(aes(yintercept=backgroudValue),linetype=5,col="red")+
      # labs(x='MNV Pattern',y='Fraction of MNVs in repetitive contexts')+
      labs(x='',y='')+
      theme_bw()+mythem2
    return(p)
  }
  
  p<-pFreq(m2_res)
  ggsave(
    filename = '../graph/one_step/2joint_freq.pdf',
    plot = p,width = 8,height = 3,
    units = 'in'
  )
}
# Repeat regions
{
  getPatternFreq<-function(df,cutoff){
    # res=list()
    df$pattern<-paste0(df$refs,'->',df$alts)
    df<-df[c(15,14,10,11)]
  # No collapsing yet
    df_matrix<-as.data.frame(table(df[1:2]),stringsAsFactors = F)
    df_matrix<-cbind(df_matrix[df_matrix$repeat.==0,c(1,3)],df_matrix[df_matrix$repeat.==1,3])
    colnames(df_matrix)<-c('pattern','notrepeat','repeat')
    df_matrix<-rbind(df_matrix,c(0,sum(df_matrix$notrepeat),sum(df_matrix$`repeat`)))
    df_matrix[nrow(df_matrix),1]<-'all'
    df_matrix$freq<-df_matrix$`repeat`/(df_matrix$`repeat`+df_matrix$notrepeat)
  # res[[1]]=df_matrix[order(-df_matrix$freq),] # optional
  # Collapse reverse-complement patterns
    df_matrix<-df_matrix[1:(nrow(df_matrix)-1),]
    df_matrix$pattern2<-apply(df_matrix[1], 1, transATCG2)
    df_matrix<-df_matrix[c(1,5,2:4)]
    df_matrix$onePattern<-apply(df_matrix[1:2],1,getPattern)
    df_matrix<-df_matrix[c(6,3,4)]
    df_list=list()
    for (i in unique(df_matrix$onePattern)) {
      df_list[[i]]<-colSums(df_matrix[df_matrix$onePattern==i,2:3])
    }
    df_matrix<-as.data.frame(do.call(rbind,df_list))
    df_matrix$pattern<-rownames(df_matrix)
    df_matrix<-df_matrix[c(3,1,2)]
    df_matrix<-rbind(df_matrix,c(0,sum(df_matrix$notrepeat),sum(df_matrix$`repeat`)))
    df_matrix[nrow(df_matrix),1]<-'all'
    df_matrix$total<-df_matrix$`repeat`+df_matrix$notrepeat
    df_matrix$freq<-df_matrix$`repeat`/df_matrix$total
    df_matrix<-df_matrix[order(-df_matrix$total),]
    flag=which(df_matrix$`repeat`==0)
    if (length(flag)>0) {
      df_matrix<-df_matrix[1:(which(df_matrix$`repeat`==0)[1]-1),]
    }
    df_matrix$realFreq<-df_matrix$freq/df_matrix[df_matrix$pattern=='all','freq']
    df_matrix<-df_matrix[df_matrix$total>cutoff,]
    return(df_matrix<-df_matrix[order(-df_matrix$freq),])
  }
  # 不考虑了连续突变的MNV 这里要跑2遍
  m2<-data[data$mnvtype==2,] 
  m2_res<-getPatternFreq(m2,100)
  m3<-data[data$mnvtype==3,] 
  m3_res<-getPatternFreq(m3,100)
  m4<-data[data$mnvtype==4 ,] 
  m4_res<-getPatternFreq(m4,100)
  m5<-data[data$mnvtype==5,] 
  m5_res<-getPatternFreq(m5,100)
  # Scatter plots
  getTop5<-function(df){
    allKey=which(df$pattern=='all')
    if (allKey>5) {
      point<-rev(df[5:1,])
    }else{
      # point<-data.frame(id=1:6,value=rev(df[1:6,'freq']),type=rev(df[1:6,'pattern']))
      point<-df[min(6,nrow(df)):1,]
      point<-point[point$pattern!='all',]
    }
    point$id<-1:nrow(point)
    # point$freq<-round(point$freq,2)*100
    return(point)
  }
  pFreq<-function(df){
    point<-getTop5(df)
    backgroudValue<-df[df$pattern=='all','freq']
    ymin=max(round(min(backgroudValue,min(point$freq)),1)-0.05,0)
    p<-ggplot(point,aes(x = id,y = freq)) +
      geom_point(aes(size=total,color=freq))+
      scale_colour_gradientn(colors = c("#FFD133", "#FF8D33", "#B0282E"))+ # color gradient
      scale_size_continuous(range=c(3,5))+ # point size range (continuous)
      # scale_size_continuous(range=c(3,5))+ # point size range (discrete)
      scale_x_continuous(breaks=c(1:nrow(point)),labels=point$pattern)+
      # scale_y_continuous(breaks=c(0:1))
      ylim(c(0,1))+
      geom_hline(aes(yintercept=backgroudValue),linetype=5,col="red")+
      # labs(x='MNV Pattern',y='Fraction of MNVs in repetitive contexts')
      labs(x='',y='')+
      theme_bw()+mythem2
    return(p)
  }
  p<-pFreq(m2_res)
  ggsave(
    filename = '../graph/repeat/2jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-pFreq(m3_res)
  ggsave(
    filename = '../graph/repeat/3jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-pFreq(m4_res)
  ggsave(
    filename = '../graph/repeat/4jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-pFreq(m5_res)
  ggsave(
    filename = '../graph/repeat/5jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  # Plot motif logos
  pMotif<-function(df_res,df){
    point<-getTop5(df_res)
    for (i in point$pattern) {
      # A,A->T,T T,T->A,A
      # i='A,A->T,T'
      ii<-transATCG2(i)
      a=str_split(i,'->',simplify = T)
      b=str_split(ii,'->',simplify = T)
      seq<-c(df[df$refs==a[1] & df$alts==a[2] & df$`repeat`==1,'ref_seq'],sapply(df[df$refs==b[1] & df$alts==b[2] & df$`repeat`==1,'ref_seq'],getTrans))
      pfm <- make_ppm(seq, ppm=FALSE)
      p<-ggseqlogo(pfm)
      ggsave(
        filename = paste0('../graph/repeat/',gsub('->','-',i),'.pdf'),
        plot = p,width = 10,height = 3,
        units = 'in'
      )
    }
  }
  pMotif(m2_res,m2)
  pMotif(m3_res,m3)
  pMotif(m4_res,m4)
  pMotif(m5_res,m5)

  # Only considering contiguous mutations below
  m2<-data[data$mnvtype==2 & data$distance==1,] 
  m2_res_continue<-getPatternFreq(m2,10)
  m3<-data[data$mnvtype==3 & data$distance==2,] 
  m3_res_continue<-getPatternFreq(m3,10)
  m4<-data[data$mnvtype==4 & data$distance==3,] 
  m4_res_continue<-getPatternFreq(m4,10)
  m5<-data[data$mnvtype==5 & data$distance==4,] 
  m5_res_continue<-getPatternFreq(m5,10)
  
  p<-pFreq(m2_res_continue)
  ggsave(
    filename = '../graph/repeat_continue/2jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-pFreq(m3_res_continue)
  ggsave(
    filename = '../graph/repeat_continue/3jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-pFreq(m4_res_continue)
  ggsave(
    filename = '../graph/repeat_continue/4jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  p<-pFreq(m5_res_continue)
  ggsave(
    filename = '../graph/repeat_continue/5jont_freq.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
  # p
  # pMotif(m5_res,m5)
  pOneMotif<-function(i,df){
    # i='A,A->T,T'
    ii<-transATCG2(i)
    a=str_split(i,'->',simplify = T)
    b=str_split(ii,'->',simplify = T)
    seq<-c(df[df$refs==a[1] & df$alts==a[2] & df$`repeat`==1,'ref_seq'],sapply(df[df$refs==b[1] & df$alts==b[2] & df$`repeat`==1,'ref_seq'],getTrans))
    pfm <- make_ppm(seq, ppm=FALSE)
    p<-ggseqlogo(pfm)
    ggsave(
      filename = paste0('../graph/repeat_continue/',gsub('->','-',i),'.pdf'),
      plot = p,width = 10,height = 3,
      units = 'in'
    )
  }
  pOneMotif('A,A->T,T',m2)
  pOneMotif('A,A,A->T,T,T',m3)
  pOneMotif('A,A,A,A->T,T,T,T',m4)
  pOneMotif('A,A,A,A,A->T,T,T,T,T',m5)
  
  pOneMotif('A,A->C,C',m2)
  pOneMotif('A,A,A->C,C,C',m3)
  pOneMotif('A,A,A,A->C,C,C,C',m4)
  pOneMotif('A,A,A,A,A->C,C,C,C,C',m5)

  pOneMotif('A,A->G,G',m2)
  pOneMotif('A,A,A->G,G,G',m3)
  pOneMotif('A,A,A,A->G,G,G,G',m4)
  pOneMotif('A,A,A,A,A->G,G,G,G,G',m5)
  
  pOneMotif('C,A->A,C',m2)
  pOneMotif('A,C->C,A',m2)
  pOneMotif('A,C,A->G,T,G',m3)
  pOneMotif('C,A,C->T,G,T',m3)
  pOneMotif('A,T,A,T,A->G,C,G,C,G',m5)

}
# ATCG frequencies by joint type
{
  getBaseNumber<-function(x){
    y<-unlist(strsplit(x,','))
    res<-c(length(which(y=='A')),length(which(y=='T')),length(which(y=='C')),length(which(y=='G')))
    return(res)
  }
  res<-apply(data['refs'], 1, getBaseNumber)
  res<-t(res)
  bases<-cbind(data$mnvtype,res) 
  bases<-as.data.frame(bases)
  colnames(bases)<-c('mnvtype','A','T','C','G')
  ATCGFreq<-list()
  for (i in 1:4) {
    ATCGFreq[[i]]<-colSums(bases[bases$mnvtype==i+1,2:5])
    ATCGFreq[[i]]<-ATCGFreq[[i]]/sum(ATCGFreq[[i]])
  }
  ATCGFreq<-do.call('rbind',ATCGFreq)
  rownames(ATCGFreq)<-paste0(2:5)
  
  ATCGFreq<-as.data.frame(ATCGFreq)
  ATCGFreq$mnvtype<-rownames(ATCGFreq)
  ATCGFreq<-reshape2::melt(ATCGFreq,'mnvtype')
  colnames(ATCGFreq)<-c('mnvtype','codon','value')
  
  p<-ggplot(ATCGFreq,aes(x=mnvtype,y=100*value, fill = codon))+
  geom_col(position="stack",width = 0.5)+ # stacked bar
    labs(x='MNV Type (N joint MNV)',y='Relative Abundance(%)')+ 
    scale_fill_manual(values = c("A" = "#2A5888", "T" = "#4396B1", "C" = "#89CEED", "G" = "#C4F5FC"))+
    theme_bw()+mythem
  ggsave(
    filename = '../graph/mnvtype_ATCG.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}
# ATCG frequencies by three mechanisms
{
  getBaseNumber<-function(x){
    y<-unlist(strsplit(x,','))
    res<-c(length(which(y=='A')),length(which(y=='T')),length(which(y=='C')),length(which(y=='G')))
    return(res)
  }
  res<-apply(data['refs'], 1, getBaseNumber)
  res<-t(res)
  bases<-cbind(data$mnvtype,res,data[12:14]) 
  bases<-as.data.frame(bases)
  colnames(bases)[1:5]<-c('mnvtype','A','T','C','G')
  
  getBasePlot<-function(bases){
    ATCGFreq<-list()
    for (i in 1:4) {
      ATCGFreq[[i]]<-colSums(bases[bases$mnvtype==i+1,2:5])
      ATCGFreq[[i]]<-ATCGFreq[[i]]/sum(ATCGFreq[[i]])
    }
    ATCGFreq<-do.call('rbind',ATCGFreq)
    rownames(ATCGFreq)<-paste0(2:5)
    
  # library(vcd)
  # spine(ATCGFreq,main='REF')
    
    library(reshape2)
    ATCGFreq<-as.data.frame(ATCGFreq)
    ATCGFreq$mnvtype<-rownames(ATCGFreq)
    ATCGFreq<-reshape2::melt(ATCGFreq,'mnvtype')
    colnames(ATCGFreq)<-c('mnvtype','codon','value')
    
    p<-ggplot(ATCGFreq,aes(x=mnvtype,y=100*value, fill = codon))+
      geom_col(position="stack",width = 0.5)+ # stack:堆叠图
      labs(x='MNV Type (N joint MNV)',y='Relative Abundance(%)')+ 
      theme_bw()+mythem
    return(p)
  }
  
  # 不同机制的不同连体的图 看看就可以了
  p1<-getBasePlot(bases[bases$snv_event==1,1:5])
  p1
  p2<-getBasePlot(bases[bases$one_step==1,1:5])
  p2  
  p3<-getBasePlot(bases[bases$`repeat`==1,1:5])
  p3
  # 绘制总的
  bases[bases$snv_event==1,1]<-1 # snv_event
  bases[bases$one_step==1,1]<-2 # one-step
  bases[bases$`repeat`==1,1]<-3  # repeat
  
  ATCGFreq<-list()
  for (i in 1:3) {
    ATCGFreq[[i]]<-colSums(bases[bases$mnvtype==i,2:5])
    ATCGFreq[[i]]<-ATCGFreq[[i]]/sum(ATCGFreq[[i]])
  }
  ATCGFreq<-do.call('rbind',ATCGFreq)
  rownames(ATCGFreq)<-c(1:3) # c('snv_event','one-step','repeat')
  
  library(reshape2)
  ATCGFreq<-as.data.frame(ATCGFreq)
  ATCGFreq$mnvtype<-rownames(ATCGFreq)
  ATCGFreq<-reshape2::melt(ATCGFreq,'mnvtype')
  colnames(ATCGFreq)<-c('mnvtype','codon','value')
  
  p<-ggplot(ATCGFreq,aes(x=mnvtype,y=100*value, fill = codon))+
  geom_col(position="stack",width = 0.5)+ # stacked bar
    labs(x='three mechanism',y='Relative Abundance(%)')+ 
    scale_fill_manual(values = c("A" = "#2A5888", "T" = "#4396B1", "C" = "#89CEED", "G" = "#C4F5FC"))+
    theme_bw()+mythem
  ggsave(
    filename = '../graph/mechanism_ATCG.pdf',
    plot = p,width = 8,height = 6,
    units = 'in'
  )
}

