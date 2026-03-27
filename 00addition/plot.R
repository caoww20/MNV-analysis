# Region distribution for each dataset
setwd('/home/caow/03mnv/analyse3/04adjustError/01region')
mythem=theme(panel.grid=element_blank(),
             plot.margin=unit(rep(2,4),'lines'),
             panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
             text = element_text(size=8,face = "bold"),
             axis.text.x=element_text(vjust=1,size=8,face = "bold", colour='black'),
             axis.title.x=element_text(vjust=1, size=8,face = "bold", color='black'),
             axis.text.y=element_text(hjust=1,size=8,face = "bold", color = 'black'),
             axis.title.y=element_text(vjust=1,size=8,face = "bold", color='black')
)
{
  fun1<-function(f){
    df<-read.table(f,header = F,stringsAsFactors = F,sep = '\t')
    colnames(df)<-c('type','num')
    return(df)
  }
  allmnv <- fun1('01num/all.txt')
  mnvRegion<-read.table('01num/region.txt',header = F,stringsAsFactors = F,sep = '\t')
  df<-cbind(fun1('01num/1000G.txt'),fun1('01num/GTEx.txt'),fun1('01num/UKB50w.txt'),fun1('01num/TCGA.txt'),fun1('01num/UKB20w.txt'))
  df<-df[c(1:2,4,6,8,10)]
  colnames(df)<-c('type','1000G','GTEx','UKB50w','TCGA','UKB20w')  
  df$id<-1:nrow(df)
  
  
  # Read SNV counts
  dfsnv<-cbind(fun1('01num/snv_1000G.region.counts.tsv'),fun1('01num/snv_GTEx.region.counts.tsv'),fun1('01num/snv_UKB50w.region.counts.tsv'),fun1('01num/snv_TCGA.region.counts.tsv'),fun1('01num/snv_UKB20w.region.counts.tsv'))
  dfsnv<-dfsnv[c(1:2,4,6,8,10)]
  colnames(dfsnv)<-c('type','1000G','GTEx','UKB50w','TCGA','UKB20w')  
  dfsnv$id<-1:nrow(dfsnv)
  
  df[,2:6] = df[,2:6]/dfsnv[,2:6]
  res<-reshape2::melt(df,id=c('type','id'))  
  
  df %>% s
  p<-ggplot(res,aes(x=id,y=value,fill=variable))+
    geom_bar(stat = 'identity',position = "dodge",show.legend = T)+
    # geom_bar(stat = 'identity',show.legend = T)+
    labs(title='',x='Annotation Type',y='Number of MNV (log10)')+
    scale_x_continuous(breaks=c(1:nrow(df)),labels=df$type)+
    theme_bw()+mythem+
    scale_fill_manual(values = c("GnomAD"="#88c6e2","1000G" = "#d02f43", "GTEx" = "#e35e42", "UKB50w" = "#eb8e68", "TCGA"="#295683","UKB20w"="#438fa9"))
  p
  
  ggsave(
    filename = './graph/bar.pdf',
    plot = p,width = 8,height = 3,
    units = 'in'
  )
  
}

# Read expected variants per region from gnomAD constraint
expected_var <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/expected_by_type.hg38.tsv')
# 1. Ensure expected_var order matches df
# left_join is safest to avoid errors from missing/extra types
df_merge <- df %>%
  left_join(expected_var, by = "type")

# 2. Define columns to include in calculations (dataset columns)
target_cols <- c("1000G", "GTEx", "UKB50w", "TCGA", "UKB20w")

# 3. Divide to compute ratios
# Store results in a new table
df_ratio <- df_merge

df_ratio[target_cols] <- sweep(df_merge[target_cols], 1, df_merge$expected_sum, "/")

# 4. Inspect results
print(head(df_ratio))
df_ratio <- df_ratio %>% select(-8)
df_ratio <- df_ratio[-15,]
res <- reshape2::melt(df_ratio,id=c('type','id')) 

ggplot(res, aes(x = id, y = value, color = variable)) +
  geom_point(size=2) +
  geom_line(lwd=1) +
  labs(title='',x='Annotation Type',y='MNV/Expected variants ratio')+
  scale_x_continuous(breaks=c(1:nrow(df_ratio)),labels=df_ratio$type)+
  theme_bw()+mythem+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_color_manual(values = c("GnomAD"="#88c6e2","1000G" = "#d02f43", "GTEx" = "#e35e42", "UKB50w" = "#eb8e68", "TCGA"="#295683","UKB20w"="#438fa9"))


# Plot five line charts
p <- ggplot(res, aes(x = id, y = value, color = variable)) +
  geom_point(size=2) +
  geom_line(lwd=1) +
  labs(title='',x='Annotation Type',y='')+
  scale_x_continuous(breaks=c(1:nrow(df)),labels=df$type)+
  theme_bw()+mythem+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_color_manual(values = c("GnomAD"="#88c6e2","1000G" = "#d02f43", "GTEx" = "#e35e42", "UKB50w" = "#eb8e68", "TCGA"="#295683","UKB20w"="#438fa9"))
ggsave(
  filename = '/home/caow/02mnv_new/mnvdensity.pdf',
  plot = p,width = 20,height = 10,
  units = 'cm'
)

library(ggplot2)
library(gridExtra)

setwd('/home/caow/03mnv/analyse3/04adjustError/01region')

# --- 1. Theme settings (unchanged) ---
mythem <- theme(
  panel.grid = element_blank(),
  # plot.margin = unit(rep(2, 4), 'lines'),
  panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
  text = element_text(size = 16, face = "bold"),
  axis.text.x = element_text(vjust = 1, size = 16, face = "bold", colour = 'black', angle = 45, hjust = 1),
  axis.title.x = element_text(vjust = 1, size = 16, face = "bold", color = 'black'),
  axis.text.y = element_text(hjust = 1, size = 16, face = "bold", color = 'black'),
  axis.title.y = element_text(vjust = 1, size = 16, face = "bold", color = 'black')
)

# --- 2. Data loading ---
fun1 <- function(f) {
  df <- read.table(f, header = F, stringsAsFactors = F, sep = '\t')
  colnames(df) <- c('type', 'num')
  return(df)
}

# Read raw MNV counts (keep original counts; do not overwrite)
df_mnv <- cbind(fun1('01num/1000G.txt'), fun1('01num/GTEx.txt'), fun1('01num/UKB50w.txt'), fun1('01num/TCGA.txt'), fun1('01num/UKB20w.txt'))
df_mnv <- df_mnv[c(1:2, 4, 6, 8, 10)]
colnames(df_mnv) <- c('type', '1000G', 'GTEx', 'UKB50w', 'TCGA', 'UKB20w')
df_mnv$id <- 1:nrow(df_mnv)

# Read raw SNV counts
df_snv <- cbind(fun1('01num/snv_1000G.region.counts.tsv'), fun1('01num/snv_GTEx.region.counts.tsv'), fun1('01num/snv_UKB50w.region.counts.tsv'), fun1('01num/snv_TCGA.region.counts.tsv'), fun1('01num/snv_UKB20w.region.counts.tsv'))
df_snv <- df_snv[c(1:2, 4, 6, 8, 10)]
colnames(df_snv) <- c('type', '1000G', 'GTEx', 'UKB50w', 'TCGA', 'UKB20w')
df_snv$id <- 1:nrow(df_snv)

# --- 3. Compute ratio data frame ---
# Copy for ratio storage
df_ratio <- df_mnv
# Compute ratio per region = MNV / SNV
df_ratio[, 2:6] <- df_mnv[, 2:6] / df_snv[, 2:6]



# --- 4. Define plotting function (includes new baseline logic) ---
plot_single_dataset <- function(data_ratio, counts_mnv, counts_snv, col_name, plot_title, line_color) {
  
  # Extract ratio data for plotting
  plot_data <- data.frame(
    id = data_ratio$id,
    type = data_ratio$type,
    value = data_ratio[[col_name]]
  )
  
  # --- Key change: compute global baseline ---
  # Definition: sum(MNV) / sum(SNV) in this dataset
  # total_mnv <- sum(counts_mnv[[col_name]], na.rm = TRUE)
  # total_snv <- sum(counts_snv[[col_name]], na.rm = TRUE)
  # baseline_value <- 0.02269008830602263
  # 
  baseline_value <- mean(plot_data$value, na.rm = TRUE)
  
  print(paste("Baseline for", col_name, ":", baseline_value)) # Print for quick check
  
  p <- ggplot(plot_data, aes(x = id, y = value)) +
    # 1. Bars (with transparency)
    geom_bar(stat = "identity", fill = line_color, alpha = 0.3, width = 0.6) +
    
    # 2. Line and points
    geom_line(color = line_color, lwd = 1) +
    geom_point(size = 2, color = line_color) +
    
    # 3. Baseline dashed line (Global Rate)
    geom_hline(yintercept = baseline_value, linetype = "dashed", color = "grey30", size = 0.8) +
    
    # 4. Axis settings
    labs(title = plot_title, y = 'MNV/SNV Ratio',x=NULL) +
    scale_x_continuous(breaks = c(1:nrow(data_ratio)), labels = gsub("_", " ", data_ratio$type)) +
    
    # 5. Theme
    theme_bw() + mythem
  
  return(p)
}

# --- 5. Run plots ---

# Plot 1000G
p_1000G <- plot_single_dataset(
  data_ratio = df_ratio,
  counts_mnv = df_mnv,   # Use raw counts for totals
  counts_snv = df_snv,   # Use raw counts for totals
  col_name = "1000G", 
  plot_title = "1000G", 
  line_color = "#d02f43"
)

# Plot GTEx
p_GTEx <- plot_single_dataset(
  data_ratio = df_ratio,
  counts_mnv = df_mnv, 
  counts_snv = df_snv,
  col_name = "GTEx", 
  plot_title = "GTEx", 
  line_color = "#e35e42"
)

# --- 6. Combine output ---
grid.arrange(p_1000G, p_GTEx, nrow = 1)

# Method B: save PDFs separately
# ggsave("plot_1000G_ratio.pdf", p_1000G, width = 5, height = 4)
# ggsave("plot_GTEx_ratio.pdf", p_GTEx, width = 5, height = 4)

mnvcategories <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/mnv_categories_1000g_gtex')
mnv1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000G_mnv_annotation')
mnvgtex <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/GTEx_mnv_annotation')
mnv1000g_anno <- mnv1000g %>% inner_join(mnvcategories,by='MNVid')
mnvgtex_anno <- mnvgtex %>% inner_join(mnvcategories,by='MNVid')

mnvmir <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/anno/mirmnv',col_names = F)
mnvmir$X2 <- 'miRNA'
colnames(mnvmir) <- colnames(mnv1000g)
mnv1000g <- rbind(mnv1000g,mnvmir)
mnvgtex <- rbind(mnvgtex,mnvmir)

mnvmir1 <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/anno/mirmnvall',col_names = F)
mnvmir1$X2 <- 'miRNA'
colnames(mnvmir1) <- colnames(mnv1000g)
mnvmir1000g <- mnvmir1 %>% inner_join(mnv1000g_pattern, by = 'MNVid') %>% select(c(1,2))
mnvmirgtex <- mnvmir1 %>% inner_join(mnvgtex_pattern, by = 'MNVid') %>% select(c(1,2))
mnv1000g <- mnv1000g %>% filter(annotation!='miRNA') %>% rbind(mnvmir1000g)
mnvgtex <- mnvgtex %>% filter(annotation!='miRNA') %>% rbind(mnvmirgtex)


setwd('/home/caow/03mnv/analyse3/05ncRNA/01pre_miRNA')
library(tidyverse)
## Read data
{
  mnv_diff<-read_tsv('02RNAFold/mnv/mnv_diff.res',col_names = F)
  snv_diff<-read_tsv('02RNAFold/snv/snv_diff.res',col_names = F)
  mnv_snv_diff<-read_tsv('02RNAFold/merge_res.txt',col_names = F)
}

# 2. Data processing (supports any number of SNV sites)
plot_data <- mnv_snv_diff %>%
  rowwise() %>% # Enable row-wise processing
  mutate(
    # Core logic:
    # 1. str_split(X4, ","): split string by comma
    # 2. unlist: flatten list to vector
    # 3. as.numeric: convert to numeric
    # 4. sum: sum values
    sum_snv_effects = sum(as.numeric(unlist(str_split(X4, ",")))),
    
    # Y axis: MNV effect
    mnv_effect = X3
  ) %>%
  ungroup() # Exit row-wise mode for performance

# 3. Define theme (unchanged)
mythem <- theme(
  panel.grid = element_blank(),
  plot.margin = unit(rep(2, 4), 'lines'),
  panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
  text = element_text(size = 8, face = "bold"),
  axis.text.x = element_text(vjust = 1, size = 8, face = "bold", colour = 'black', hjust = 1),
  axis.title.x = element_text(vjust = 1, size = 8, face = "bold", color = 'black'),
  axis.text.y = element_text(hjust = 1, size = 8, face = "bold", color = 'black'),
  axis.title.y = element_text(vjust = 1, size = 8, face = "bold", color = 'black')
)
library(ggpubr)
# 4. Plot
p <- ggplot(plot_data, aes(x = sum_snv_effects, y = mnv_effect)) +
  geom_point(color = "#377EB8", alpha = 0.7, size = 1.5) + 
  
  # Add diagonal (y=x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#E41A1C", size = 0.8) +
  
  # Axis labels
  # Note: use a more general label (not specific to SNV1+SNV2)
  labs(
    x = expression(bold(Sum~of~individual~SNV~effects~(Sigma*Delta*MFE[SNV]))),
    y = expression(bold(MNV~effect~(Delta*MFE[MNV])))
  ) +
  
  theme_bw() + 
  mythem
p_with_cor <- p + 
  stat_cor(
    method = "pearson",      # Use Pearson correlation
    label.x.npc = 0.05,      # Horizontal label position (0-1, 0.05 is left)
    label.y.npc = 0.95,      # Vertical label position (0-1, 0.95 is top)
    # Or use absolute coordinates: label.x = -10, label.y = 10
    size = 4,                # Font size
    fontface = "bold"        # Bold text
  )
# 5. Display and save
# ggsave("MNV_vs_MultiSNVs_scatter.pdf", p, width = 5, height = 5)

threshold_value <- 1.0 

# 2. Compute deviations and proportions
stats_result <- plot_data %>%
  mutate(
    # Residuals: observed - expected
    deviation = mnv_effect - sum_snv_effects,
    
    # Absolute value for thresholding
    # abs() ignores direction; only magnitude matters
    is_non_additive = abs(deviation) > threshold_value
  )

# 3. Summarize counts
total_count <- nrow(stats_result)
non_additive_count <- sum(stats_result$is_non_additive)
percentage <- round((non_additive_count / total_count) * 100, 2)

mnvcategories <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/mnv_categories_1000g_gtex')
mnv1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000G_mnv_annotation')
mnvgtex <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/GTEx_mnv_annotation')
mnv1000g_anno <- mnv1000g %>% inner_join(mnvcategories,by='MNVid')
mnvgtex_anno <- mnvgtex %>% inner_join(mnvcategories,by='MNVid')

custom_colors <- c(
  "snv_event" = "#A83232",   # Deep brick red
  "one_step"  = "#E8A03D",   # Orange
  "repeat"    = "#2C4B8E"    # Deep blue
)

# B. Define label text (for figure 2)
custom_labels <- c(
  "snv_event" = "SNV event",
  "one_step"  = "One-step",
  "repeat"    = "Slippage at Repeat"
)

my_order <- c('UTR5','exon','splice','intron','UTR3','lncRNA','miRNA','circRNA','piRNA','ATAC','conserved_element','enhancer','miRBS','TFBS','intergenic')
# 2. Data cleaning and transformation (core step)
plot_data <- mnv1000g_anno %>%
  # Step 1: split semicolon-separated annotations into multiple rows
  # Example: "lncRNA;piRNA" becomes two rows: lncRNA and piRNA
  separate_rows(annotation, sep = ";") %>%
  
  # Step 2: pivot category columns to long format
  # Keep only rows with value 1 (the MNV belongs to that category)
  pivot_longer(
    cols = c(snv_event, one_step, `repeat`), 
    names_to = "category", 
    values_to = "is_present"
  ) %>%
  
  # Filter out rows with value 0 (not in that category)
  filter(is_present == 1)

plot_data <- plot_data %>%
  mutate(annotation = factor(annotation, levels = my_order))

my_order <- c('UTR5','exon','splice','intron','UTR3','lncRNA','miRNA','circRNA','piRNA','ATAC','conserved_element','enhancer','miRBS','TFBS','intergenic','TSS')
# 2. Data cleaning and transformation (core step)
plot_data <- mnv1000g_anno %>%
  # Step 1: split semicolon-separated annotations into multiple rows
  # Example: "lncRNA;piRNA" becomes two rows: lncRNA and piRNA
  separate_rows(annotation, sep = ";") %>%
  
  # Step 2: pivot category columns to long format
  # Keep only rows with value 1 (the MNV belongs to that category)
  pivot_longer(
    cols = c(snv_event, one_step, `repeat`), 
    names_to = "category", 
    values_to = "is_present"
  ) %>%
  
  # Filter out rows with value 0 (not in that category)
  filter(is_present == 1)

# Compare proportions of MNV categories across functional regions
plot_data <- plot_data %>%
  mutate(annotation = factor(annotation, levels = my_order))

mnvtypenum <- plot_data %>% group_by(annotation,category) %>% summarise(n=n())
colnames(mnvtypenum) <- c('type','category','1000G')

plot_data1_gtex <- mnvgtex_anno %>% distinct() %>% 
  # Step 1: split semicolon-separated annotations into multiple rows
  # Example: "lncRNA;piRNA" becomes two rows: lncRNA and piRNA
  separate_rows(annotation, sep = ";") %>%
  
  # Step 2: pivot category columns to long format
  # Keep only rows with value 1 (the MNV belongs to that category)
  pivot_longer(
    cols = c(snv_event, one_step, `repeat`), 
    names_to = "category", 
    values_to = "is_present"
  ) %>%
  
  # Filter out rows with value 0 (not in that category)
  filter(is_present == 1)

plot_data1_gtex <- plot_data1_gtex %>%
  mutate(annotation = factor(annotation, levels = my_order))

mnvtypenum_gtex <- plot_data1_gtex %>% group_by(annotation,category) %>% summarise(n=n())
colnames(mnvtypenum_gtex) <- c('type','category','GTEx')

mnvtypenum <-  mnvtypenum %>% left_join(mnvtypenum_gtex,by=c('type','category'))
mnvtypenum1 <-  mnvtypenum %>% filter(category=='one_step') %>% select(-2)
mnvtypenum2 <-  mnvtypenum %>% filter(category=='repeat') %>% select(-2)
mnvtypenum3 <-  mnvtypenum %>% filter(category=='snv_event') %>% select(-2)


# Read raw SNV counts
df_snv <- cbind(fun1('01num/snv_1000G.region.counts.tsv'), fun1('01num/snv_GTEx.region.counts.tsv'), fun1('01num/snv_UKB50w.region.counts.tsv'), fun1('01num/snv_TCGA.region.counts.tsv'), fun1('01num/snv_UKB20w.region.counts.tsv'))
df_snv <- df_snv[c(1:2, 4, 6, 8, 10)]
colnames(df_snv) <- c('type', '1000G', 'GTEx', 'UKB50w', 'TCGA', 'UKB20w')
df_snv$id <- 1:nrow(df_snv)

# --- 3. Compute ratio data frame ---
# Copy for ratio storage
# df_mnv <- df_mnv[-17,]

df_ratio1 <- df_mnv
df_ratio2 <- df_mnv
df_ratio3 <- df_mnv

# df_snv <- df_snv[-16,]

# Compute ratio per region = MNV / SNV
df_ratio1[, 2:3] <- mnvtypenum1[-17, 2:3] / df_snv[, 2:3]
df_ratio2[, 2:3] <- mnvtypenum2[-17, 2:3] / df_snv[, 2:3]
df_ratio3[, 2:3] <- mnvtypenum3[-17, 2:3] / df_snv[, 2:3]

plot_single_dataset1 <- function(data_ratio, counts_mnv, counts_snv, col_name, plot_title, line_color, baseline) {
  
  # Extract ratio data for plotting
  plot_data <- data.frame(
    id = data_ratio$id,
    type = data_ratio$type,
    value = data_ratio[[col_name]]
  )
  
  # --- Key change: compute global baseline ---
  # Definition: sum(MNV) / sum(SNV) in this dataset
  # total_mnv <- sum(counts_mnv[[col_name]], na.rm = TRUE)
  # total_snv <- sum(counts_snv[[col_name]], na.rm = TRUE)
  # baseline_value <- 0.02269008830602263
  # 
  baseline_value <- baseline
  
  print(paste("Baseline for", col_name, ":", baseline_value)) # Print for quick check
  
  p <- ggplot(plot_data, aes(x = id, y = value)) +
    # 1. Bars (with transparency)
    geom_bar(stat = "identity", fill = line_color, alpha = 0.3, width = 0.6) +
    
    # 2. Line and points
    geom_line(color = line_color, lwd = 1) +
    geom_point(size = 2, color = line_color) +
    
    # 3. Baseline dashed line (Global Rate)
    geom_hline(yintercept = baseline_value, linetype = "dashed", color = "grey30", size = 0.8) +
    
    # 4. Axis settings
    labs(title = plot_title, x = 'Genomic Region', y = 'MNV-to-SNV Ratio') +
    scale_x_continuous(breaks = c(1:nrow(data_ratio)), labels = data_ratio$type) +
    
    # 5. Theme
    theme_bw() + mythem
  
  return(p)
}

# Plot 1000G
p_1000G1 <- plot_single_dataset1(
  data_ratio = df_ratio1,
  counts_mnv = mnvtypenum1,   # Use raw counts for totals
  counts_snv = df_snv,   # Use raw counts for totals
  col_name = "1000G", 
  plot_title = "1000G Dataset", 
  line_color = "#d02f43",baseline = 0.041212660155921
)

# Plot GTEx
p_GTEx1 <- plot_single_dataset1(
  data_ratio = df_ratio1,
  counts_mnv = mnvtypenum1, 
  counts_snv = df_snv,
  col_name = "GTEx", 
  plot_title = "GTEx Dataset", 
  line_color = "#e35e42",baseline = 0.027862315966804
)

p_1000G2 <- plot_single_dataset1(
  data_ratio = df_ratio2,
  counts_mnv = mnvtypenum1,   # Use raw counts for totals
  counts_snv = df_snv,   # Use raw counts for totals
  col_name = "1000G", 
  plot_title = "1000G Dataset", 
  line_color = "#d02f43",baseline = 0.041212660155921
)

# Plot GTEx
p_GTEx2 <- plot_single_dataset1(
  data_ratio = df_ratio2,
  counts_mnv = mnvtypenum1, 
  counts_snv = df_snv,
  col_name = "GTEx", 
  plot_title = "GTEx Dataset", 
  line_color = "#e35e42",baseline = 0.027862315966804
)

p_1000G3 <- plot_single_dataset1(
  data_ratio = df_ratio3,
  counts_mnv = mnvtypenum1,   # Use raw counts for totals
  counts_snv = df_snv,   # Use raw counts for totals
  col_name = "1000G", 
  plot_title = "1000G Dataset", 
  line_color = "#d02f43",baseline = 0.041212660155921
)

# Plot GTEx
p_GTEx3 <- plot_single_dataset1(
  data_ratio = df_ratio3,
  counts_mnv = mnvtypenum1, 
  counts_snv = df_snv,
  col_name = "GTEx", 
  plot_title = "GTEx Dataset", 
  line_color = "#e35e42",baseline = 0.027862315966804
)
library(grid)
library(gridExtra)

# Define title style for each row
title_style <- function(label) {
  textGrob(label, gp = gpar(fontsize = 12, fontface = "bold"))
}

# Step 1: combine each row with its mechanism label
# Row 1: SNV event
row1 <- arrangeGrob(p_1000G1, p_GTEx1, ncol = 2, 
                    top = title_style("pol-zeta-mediated one-step"))

# Row 2: pol-zeta-mediated one-step
row2 <- arrangeGrob(p_1000G2, p_GTEx2, ncol = 2, 
                    top = title_style("repeat slippage"))

# Row 3: repeat slippage
row3 <- arrangeGrob(p_1000G3, p_GTEx3, ncol = 2, 
                    top = title_style("SNV event"))

# Step 2: stack three rows vertically
grid.arrange(row1, row2, row3, nrow = 3)


prepare_plot_data <- function(df1, df2, df3, dataset_col) {
  # Extract Pol-zeta
  d1 <- df1[, c("id", "type", dataset_col)]
  colnames(d1)[3] <- "value"
  d1$Type <- "Pol-zeta-mediated one-step"
  
  # Extract Repeat slippage
  d2 <- df2[, c("id", "type", dataset_col)]
  colnames(d2)[3] <- "value"
  d2$Type <- "Repeat slippage"
  
  # Extract SNV event
  d3 <- df3[, c("id", "type", dataset_col)]
  colnames(d3)[3] <- "value"
  d3$Type <- "SNV event"
  
  # Merge data
  plot_data <- rbind(d1, d2, d3)
  
  # Ensure X-axis order: factor type ordered by id
  # Assume id is 1, 2, 3... matching original row order
  order_levels <- unique(d1$type[order(d1$id)])
  plot_data$type <- factor(plot_data$type, levels = order_levels)
  
  # Ensure legend order for Mechanism
  plot_data$Type <- factor(plot_data$Type, 
                                levels = c("SNV event", "Pol-zeta-mediated one-step", "Repeat slippage"))
  
  return(plot_data)
}

# =======================================================
# 2. Define plotting function (multi-mechanism)
# =======================================================
plot_lines_final <- function(plot_data, plot_title, baseline_value, show_legend = TRUE, show_x_title = TRUE) {
  
  print(paste("Drawing:", plot_title, "| Baseline:", baseline_value)) 
  
  # Base plot object
  p <- ggplot(plot_data, aes(x = type, y = value, color = Type, group = Type)) +
    
    # 1. Baseline
    # geom_hline(yintercept = baseline_value, linetype = "dashed", color = "grey40", size = 0.8) +
    
    # 2. Line and points
    geom_line(lwd=1, alpha = 0.9) + 
    geom_point(size = 2) +

    # 3. Color settings
    scale_color_manual(values = c("SNV event" = "#A83232", 
                                  "Pol-zeta-mediated one-step" = "#E8A03D", 
                                  "Repeat slippage" = "#2C4B8E")) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    # 4. Title and axis label control
    labs(title = plot_title, 
         # Show or hide X-axis title based on parameter
         x = if(show_x_title) 'Genomic Region' else NULL, 
         y = 'MNV/SNV Ratio') +
    

    # 5. Theme settings
    theme_bw() + mythem+theme(legend.position = if(show_legend) "top" else "none")
    # theme(
    #   # --- Key change: center title ---
    #   plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    #   
    #   # Axis labels
    #   axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    #   axis.title.y = element_text(size = 12, face = "bold"),
    #   
    #   # --- Key change: legend position control ---
    #   # If show_legend is TRUE, place at top; otherwise hide
    #   legend.position = if(show_legend) "top" else "none",
    #   legend.title = element_blank(),
    #   legend.text = element_text(size = 11, face = "bold"),
    #   
    #   panel.grid.minor = element_blank()
    #   # + mythem (if needed)
    # )
  
  return(p)
}

# =======================================================
# 3. Run plots
# =======================================================

# --- Prepare 1000G data ---
# Note: assumes df_ratio1, df_ratio2, df_ratio3 are already computed
data_1000G <- prepare_plot_data(df_ratio1, df_ratio2, df_ratio3, "1000G")

# --- Prepare GTEx data ---
data_GTEx <- prepare_plot_data(df_ratio1, df_ratio2, df_ratio3, "GTEx")

library(ggpubr) # install.packages("ggpubr")

# 1. Ensure show_legend = TRUE in both plots
# (ggpubr will extract one legend and hide the others)
p_1000G_with_leg <- plot_lines_final(data_1000G, "1000G", 0.0412, show_legend = TRUE, show_x_title = FALSE)
p_GTEx_with_leg  <- plot_lines_final(data_GTEx,  "GTEx",  0.0278, show_legend = TRUE, show_x_title = FALSE)

# 2. Arrange with shared legend
p2 <- ggarrange(
  p_1000G_with_leg, 
  p_GTEx_with_leg, 
  ncol = 2, 
  nrow = 1,
  common.legend = TRUE, # Enable shared legend
  legend = "bottom"     # Place legend at bottom
)
p1 <- grid.arrange(p_1000G, p_GTEx, nrow = 1)
grid.arrange(p1, p2, ncol = 1)

library(cowplot)
library(ggpubr)

# --- Step 1: build lower row with shared legend (Row 2) ---
# Use ggarrange to extract shared legend
p2 <- ggarrange(
  p_1000G_with_leg, 
  p_GTEx_with_leg, 
  ncol = 2, 
  nrow = 1,
  common.legend = TRUE, 
  legend = "bottom"
)

# --- Step 2: build upper row (Row 1) ---
# plot_grid is preferred over grid.arrange for a ggplot object
p1 <- plot_grid(
  p_1000G, 
  p_GTEx, 
  ncol = 2, 
  align = "h" # Horizontal alignment
)

# --- Step 3: final composition and spacing ---
final_plot <- plot_grid(
  p1, p2, 
  ncol = 1, 
  
  # 1. Key: align axes between rows
  align = "v", 
  axis = "lr",
  
  # 2. Key: adjust row height ratio
  # Lower row has legend and X-axis title, so give it more space (e.g. 1:1.2)
  rel_heights = c(1, 1.1),
  
  # 3. Key: increase spacing between subplots (optional)
  # scale < 1 slightly shrinks plots to add whitespace
  scale = c(1.07, 1.07) 
)

final_plot

# --- Plot ---
# Baseline can be passed in or hard-coded in the function
# 1000G plot (Baseline: 0.0412...)
p_top <- plot_lines_final(
  plot_data = data_1000G, 
  plot_title = "1000G Dataset", 
  baseline_value = 0.041212660155921,
  show_legend = TRUE,     # Show legend
  show_x_title = TRUE    # Hide X-axis title
)

# --- Lower plot (GTEx) ---
# Hide legend (show_legend = FALSE) since it is already shown above
# Keep X-axis title (show_x_title = TRUE)
p_bottom <- plot_lines_final(
  plot_data = data_GTEx, 
  plot_title = "GTEx Dataset", 
  baseline_value = 0.027862315966804,
  show_legend = FALSE,    # Hide legend
  show_x_title = TRUE     # Show X-axis title
)

grid.arrange(p_top, p_bottom, nrow = 1)


library(cowplot) # Recommended for arranging plots

plot_grid(p_top, p_bottom, 
          ncol = 1, 
          align = "h", # Align axes vertically
          rel_heights = c(1.2, 1.1))


# 3. Plot (percentage stacked bar)
ggplot(plot_data, aes(x = annotation, fill = category)) +
  # position = "fill" computes proportions and stacks
  geom_bar(position = "stack", width = 0.7) + 
  
  # Add percentage labels (optional for clarity)
  geom_text(stat = "count", 
            aes(label = scales::percent(after_stat(count) / tapply(after_stat(count), after_stat(x), sum)[after_stat(x)], accuracy = 1)),
            position = position_fill(vjust = 0.5), 
            size = 3, color = "black") +
  
  # Axis and color adjustments
  scale_y_continuous(labels = scales::percent) + # Y axis in percent
  scale_fill_manual(
    values = custom_colors, # Use defined color palette
    labels = custom_labels, # Use defined labels
    name = NULL             # Remove legend title (or set name = "Variant Type")
  ) +
  labs(
    title = "Proportion of MNV Categories by Annotation Region",
    x = "Genomic Region (Annotation)",
    y = "Proportion",
    fill = "MNV Type"
  ) +
  theme_bw() + mythem+
  # Rotate long x-axis labels if needed
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# To show counts instead of proportions, change position = "fill" to position = "stack"


mnv1000g_pattern <- read_tsv('~/03mnv/analyse3/00addition/data/mnv1000g_pattern',col_names = F)
colnames(mnv1000g_pattern) <- c('MNVid','ref','alt','dis')

# Create MNV classification data frame
mnv_data <- read.csv(text = "
refs,alts,class
AA,CC,Tv
AA,CG,TiTv
AA,CT,Tv
AA,GC,TiTv
AA,GG,Ti
AA,GT,TiTv
AA,TC,Tv
AA,TG,TiTv
AA,TT,repeat
AC,CA,Tv
AC,CG,Tv
AC,CT,CpG_TiTv
AC,GA,TiTv
AC,GG,TiTv
AC,GT,CpG_Ti
AC,TA,Tv
AC,TG,Tv
AC,TT,CpG_TiTv
AG,CA,CpG_TiTv
AG,CC,Tv
AG,CT,Tv
AG,GA,Ti
AG,GC,TiTv
AG,GT,TiTv
AG,TA,TiTv
AG,TC,Tv
AG,TT,Tv
AT,CA,Tv
AT,CC,TiTv
AT,CG,Tv
AT,GA,TiTv
AT,GC,Ti
AT,TA,repeat
CA,AC,Tv
CA,AG,TiTv
CA,AT,Tv
CA,GC,Tv
CA,GG,TiTv
CA,GT,Tv
CA,TC,TiTv
CA,TG,CpG_Ti
CA,TT,TiTv
CC,AA,Tv
CC,AG,Tv
CC,AT,CpG_TiTv
CC,GA,Tv
CC,GG,Tv
CC,GT,CpG_TiTv
CC,TA,TiTv
CC,TG,CpG_TiTv
CC,TT,CpG_Ti
CG,AA,CpG_TiTv
CG,AC,Tv
CG,AT,Tv
CG,GA,CpG_TiTv
CG,GC,Tv
CG,TA,CpG_Ti
GA,AC,CpG_TiTv
GA,AG,CpG_Ti
GA,AT,CpG_TiTv
GA,CC,Tv
GA,CG,TiTv
GA,CT,Tv
GA,TC,Tv
GA,TG,TiTv
GA,TT,pol zeta
GC,AA,pol zeta
GC,AG,CpG_TiTv
GC,AT,CpG_Ti
GC,CA,Tv
GC,CG,Tv
GC,TA,Tv
TA,AC,Tv
TA,AG,TiTv
TA,AT,repeat
TA,CC,TiTv
TA,CG,Ti
TA,GC,Tv
", strip.white = TRUE, stringsAsFactors = FALSE)

# Add mutation column (format like AA->CC)
mnv_data$mutation <- paste(mnv_data$refs, "->", mnv_data$alts, sep="")

# Reorder columns to match the figure layout
mnv_data <- mnv_data[, c("mutation", "refs", "alts", "class")]

# --- build ref_dict ---
basemap <- c(A = 1, G = 3, C = 2, T = 4)
elements <- c("A", "C", "G", "T")

combn_mat <- combn(elements, 2)  # 2 x 6
ref_dict <- character(0)
for (idx in seq_len(ncol(combn_mat))) {
  a <- combn_mat[, idx]
  ref_dict <- c(ref_dict,
                paste(a, collapse = ","),
                paste(rev(a), collapse = ","))
}

# --- build all variant keys like 'A,T:G,A' ---
variant_keys <- character(0)
for (i in ref_dict) {
  for (j in ref_dict) {
    # i like "A,C" -> i[1] is "A", i[3] is "C" in R substr sense
    t <- paste0(substr(i, 1, 1), ",", substr(j, 1, 1), ":",
                substr(i, 3, 3), ",", substr(j, 3, 3))
    variant_keys <- c(variant_keys, t)
  }
}
variant_keys <- unique(variant_keys)

# (Optional) If you also want r / r1 / r_all structures, build lists here:
# r  <- setNames(replicate(length(variant_keys), integer(3203), simplify = FALSE), variant_keys)
# r1 <- setNames(replicate(length(variant_keys), integer(3203), simplify = FALSE), variant_keys)

# --- reverse complement function ---
reverse_complement_variant <- function(variant) {
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  parts <- strsplit(variant, ":", fixed = TRUE)[[1]]
  ref <- strsplit(parts[1], ",", fixed = TRUE)[[1]]
  alt <- strsplit(parts[2], ",", fixed = TRUE)[[1]]
  
  ref_comp <- unname(complement[ref])
  alt_comp <- unname(complement[alt])
  
  paste0(ref_comp[2], ",", ref_comp[1], ":", alt_comp[2], ",", alt_comp[1])
}

# --- build mapdict ---
mapdict <- character(0)

for (i in variant_keys) {
  rc <- reverse_complement_variant(i)
  
  if (!identical(i, rc)) {
    # Compare i and rc bases at positions 1,3,5,7 (like Python range(0,8,2))
    decided <- FALSE
    for (pos in c(1, 3, 5, 7)) {
      bi <- substr(i, pos, pos)
      bj <- substr(rc, pos, pos)
      
      if (basemap[bi] < basemap[bj]) {
        mapdict[i] <- rc
        decided <- TRUE
        break
      } else if (basemap[bi] > basemap[bj]) {
        mapdict[rc] <- i
        decided <- TRUE
        break
      } else {
        next
      }
    }
    # decided should always be set; if fully equal then i==rc (already excluded)
  } else {
    mapdict[i] <- i
  }
}

# mapdict is the final result
# Example: view first mappings
head(mapdict)

mnv1000g_2joint_pattern <- mnv1000g %>% inner_join(mnv1000g_pattern,by='MNVid') %>% filter(nchar(ref)==3&dis==1)

# --- Step 1: build full symmetric lookup table ---
# mapdict is {smaller: larger} or {self: self}
# We need a table mapping larger -> smaller, and smaller -> itself
smaller_keys <- names(mapdict)
larger_values <- unname(mapdict)

# Build a full lookup vector: map all raw patterns to representative patterns
full_lookup_vec <- setNames(smaller_keys, smaller_keys) # smaller -> smaller
full_lookup_vec <- c(full_lookup_vec, setNames(smaller_keys, larger_values)) # larger -> smaller
# De-duplicate (avoid adding self-complement patterns twice)
full_lookup_vec <- full_lookup_vec[!duplicated(names(full_lookup_vec))]

# Add promoter and TSS to mnv1000g
promotermnv1000G <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/anno/promotermnv1000G',col_names = F) %>% 
  mutate(X2='Promoter')
promotermnvGTEx <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/anno/promotermnvGTEx',col_names = F) %>% 
  mutate(X2='Promoter')
tss1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000gMNV_to_TSS.tsv',col_names = F) %>% 
  mutate(X2='TSS')
tssGTEx <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/GTExMNV_to_TSS.tsv',col_names = F) %>% 
  mutate(X2='TSS')
colnames(promotermnv1000G) <- colnames(mnv1000g)
colnames(promotermnvGTEx) <- colnames(mnv1000g)
colnames(tss1000g) <- colnames(mnv1000g)
colnames(tssGTEx) <- colnames(mnv1000g)
mnv1000g <- rbind(mnv1000g,promotermnv1000G,tss1000g)
mnvgtex <- rbind(mnvgtex,promotermnvGTEx,tssGTEx)

# --- Step 2: process mnv1000g data ---
mnv1000g_classified <- mnv1000g %>% 
  inner_join(mnv1000g_pattern, by = 'MNVid') %>% 
  # Filter 2-joint and adjacent MNVs
  filter(nchar(ref) == 3 & dis == 1) %>% 
  mutate(
    # 1. Build raw mapdict format, e.g. "C,G:T,A"
    raw_key = paste0(ref, ":", alt)
  ) %>% 
  mutate(
    # 2. Use lookup table to get representative pattern
    # This is fast via vectorized replacement
    rep_key = full_lookup_vec[raw_key]
  ) %>% 
  mutate(
    # 3. Convert "A,A:C,C" to "AA->CC" to match mnv_data
    match_pattern = str_replace_all(rep_key, ",", "") %>% str_replace(":", "->")
  ) %>% 
  # 4. Join mnv_data to get class
  left_join(mnv_data, by = c("match_pattern" = "mutation"))

mnvgtex_pattern <- read_tsv('~/03mnv/analyse3/00addition/data/mnvGTEx_pattern',col_names = F)
colnames(mnvgtex_pattern) <- c('MNVid','ref','alt','dis')
mnvgtex_2joint_pattern <- mnvgtex %>% inner_join(mnvgtex_pattern,by='MNVid') %>% filter(nchar(ref)==3&dis==1)

mnvgtex_classified <- mnvgtex %>% 
  inner_join(mnvgtex_pattern, by = 'MNVid') %>% 
  # Filter 2-joint and adjacent MNVs
  filter(nchar(ref) == 3 & dis == 1) %>% 
  mutate(
    # 1. Build raw mapdict format, e.g. "C,G:T,A"
    raw_key = paste0(ref, ":", alt)
  ) %>% 
  mutate(
    # 2. Use lookup table to get representative pattern
    # This is fast via vectorized replacement
    rep_key = full_lookup_vec[raw_key]
  ) %>% 
  mutate(
    # 3. Convert "A,A:C,C" to "AA->CC" to match mnv_data
    match_pattern = str_replace_all(rep_key, ",", "") %>% str_replace(":", "->")
  ) %>% 
  # 4. Join mnv_data to get class
  left_join(mnv_data, by = c("match_pattern" = "mutation"))

# --- Step 3: review summary ---
# Count per class
class_summary <- mnv1000g_classified %>%
  group_by(class) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

print(class_summary)

# 1. Define mutation pattern lists (per requirement)
cpg_list     <- c("GA->AG", "CC->TT", "AC->GT", "CA->TG")
polzeta_list <- c("GA->TT", "GC->AA")
rep_list     <- c("AC->CA", "CA->AC", "AT->TA", "TA->AT", "AA->TT", "AA->CC")
tv_list      <- c("TA->CG", "CG->AT", "AT->CG", "CG->GC", "GC->CG", "CG->AC")

# 2. Define region order
region_order <- c("enhancer", "Promoter", "TSS", "UTR5", "exon", "intron", "UTR3")

# 3. Data processing
plot_df <- mnv1000g_classified %>%
  # Handle multi-annotations: split by semicolon
  separate_rows(annotation, sep = ";") %>%
  
  # --- Key change: map mechanism by specific pattern ---
  # Note: uses final/match pattern in "AA->CC" format
  mutate(mechanism = case_when(
    match_pattern %in% cpg_list     ~ "Ti at CpG",
    match_pattern %in% polzeta_list ~ "Pol-zeta error",
    match_pattern %in% rep_list     ~ "Slippage at repeat",
    match_pattern %in% tv_list      ~ "Tv combination"
  )) %>%
  
  # Filter target regions
  filter(annotation %in% region_order) %>%
  
  # Count fractions by region and mechanism
  group_by(annotation, mechanism) %>%
  summarise(count = n(), .groups = 'drop') %>%
  
  group_by(annotation) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup() %>%
  
  # Set factor order
  mutate(annotation = factor(annotation, levels = region_order)) %>% filter(!is.na(mechanism))

# 4. Plot
ggplot(plot_df, aes(x = annotation, y = fraction, color = mechanism, group = mechanism)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +
  
  # Color scheme (match figure)
  scale_color_manual(values = c(
    "Ti at CpG"          = "#56B4E9", # Light blue
    "Pol-zeta error"     = "#9970AB", # Purple
    "Slippage at repeat" = "#8C510A", # Brown
    "Tv combination"     = "#D6604D" # Red
  )) +
  
  labs(
    x = NULL, 
    y = "Estimated fraction of origin",
    color = NULL
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "grey80"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.5)
  )

plot_df1 <- mnvgtex_classified %>%
  # Handle multi-annotations: split by semicolon
  separate_rows(annotation, sep = ";") %>%
  
  # --- Key change: map mechanism by specific pattern ---
  # Note: uses final/match pattern in "AA->CC" format
  mutate(mechanism = case_when(
    match_pattern %in% cpg_list     ~ "Ti at CpG",
    match_pattern %in% polzeta_list ~ "Pol-zeta error",
    match_pattern %in% rep_list     ~ "Slippage at repeat",
    match_pattern %in% tv_list      ~ "Tv combination"
  )) %>%
  
  # Filter target regions
  filter(annotation %in% region_order) %>%
  
  # Count fractions by region and mechanism
  group_by(annotation, mechanism) %>%
  summarise(count = n(), .groups = 'drop') %>%
  
  group_by(annotation) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup() %>%
  
  # Set factor order
  mutate(annotation = factor(annotation, levels = region_order)) %>% filter(!is.na(mechanism))

# 4. Plot
ggplot(plot_df1, aes(x = annotation, y = fraction, color = mechanism, group = mechanism)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +
  
  # Color scheme (match figure)
  scale_color_manual(values = c(
    "Ti at CpG"          = "#56B4E9", # Light blue
    "Pol-zeta error"     = "#9970AB", # Purple
    "Slippage at repeat" = "#8C510A", # Brown
    "Tv combination"     = "#D6604D" # Red
  )) +
  
  labs(
    x = NULL, 
    y = "Estimated fraction of origin",
    color = NULL
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "grey80"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.5)
  )

mnvmerge_pattern <- read_tsv('~/03mnv/analyse3/00addition/data/mnv1000gGTEx_pattern')
colnames(mnvmerge_pattern) <- c('MNVid','ref','alt','dis')
mnv1000ggtex <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000G_GTEx_annotation.tsv')


mnv1000ggtex <- rbind(mnv1000ggtex,mnvmir)

promotermnvall <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/promotermnv_all',col_names = F) %>% 
  mutate(X2='Promoter')
colnames(promotermnvall) <- colnames(mnv1000g)

tssall <- rbind(tss1000g,tssGTEx) %>% distinct()
mnv1000ggtex<- rbind(mnv1000ggtex,promotermnvall,tssall)

utr5mnv <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/mnv_in_utr5_ids.txt',col_names = F) %>% 
  mutate(X2='UTR5')
colnames(utr5mnv) <- colnames(mnv1000g)
mnv1000ggtex <-  rbind(mnv1000ggtex,utr5mnv)

mnvmerge_classified <- mnv1000ggtex %>% 
  inner_join(mnvmerge_pattern, by = 'MNVid') %>% 
  # Filter 2-joint and adjacent MNVs
  filter(nchar(ref) == 3 & dis == 1) %>% 
  mutate(
    # 1. Build raw mapdict format, e.g. "C,G:T,A"
    raw_key = paste0(ref, ":", alt)
  ) %>% 
  mutate(
    # 2. Use lookup table to get representative pattern
    # This is fast via vectorized replacement
    rep_key = full_lookup_vec[raw_key]
  ) %>% 
  mutate(
    # 3. Convert "A,A:C,C" to "AA->CC" to match mnv_data
    match_pattern = str_replace_all(rep_key, ",", "") %>% str_replace(":", "->")
  ) %>% 
  # 4. Join mnv_data to get class
  left_join(mnv_data, by = c("match_pattern" = "mutation"))
plot_df2 <- mnvmerge_classified %>%
  # Handle multi-annotations: split by semicolon
  separate_rows(annotation, sep = ";") %>% distinct() %>% 
  
  # --- Key change: map mechanism by specific pattern ---
  # Note: uses final/match pattern in "AA->CC" format
  mutate(mechanism = case_when(
    match_pattern %in% cpg_list     ~ "Ti at CpG",
    match_pattern %in% polzeta_list ~ "Pol-zeta error",
    match_pattern %in% rep_list     ~ "Slippage at repeat",
    match_pattern %in% tv_list      ~ "Tv combination"
  )) %>%
  
  # Filter target regions
  filter(annotation %in% region_order) %>%
  
  # Count fractions by region and mechanism
  group_by(annotation, mechanism) %>%
  summarise(count = n(), .groups = 'drop') %>%
  
  group_by(annotation) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup() %>%
  
  # Set factor order
  mutate(annotation = factor(annotation, levels = region_order)) %>% filter(!is.na(mechanism))

# 4. Plot
ggplot(plot_df2, aes(x = annotation, y = fraction, color = mechanism, group = mechanism)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +
  
  # Color scheme (match figure)
  scale_color_manual(values = c(
    "Ti at CpG"          = "#56B4E9", # Light blue
    "Pol-zeta error"     = "#9970AB", # Purple
    "Slippage at repeat" = "#8C510A", # Brown
    "Tv combination"     = "#D6604D" # Red
  )) +
  
  labs(
    x = NULL, 
    y = "Estimated fraction of origin",
    color = NULL
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "grey80"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.5)
  )



# --- 1. Prepare plot data ---
# Assume mnv1000g_classified already has 'mechanism' and one-step logic
# We need the one-step fraction per pattern

plot_pattern_df <- mnv1000g_classified %>%
  # Determine one-step (example logic; compute from AC columns in your data)
  # mutate(is_onestep = ifelse(AC_mnv / AC_snv_min > 0.9 & (AC_snv_max - AC_snv_min)/AC_snv_min < 0.1, TRUE, FALSE)) %>%
  
  group_by(mutation_pattern, mechanism) %>%
  summarise(
    onestep_fraction = mean(is_onestep, na.rm = TRUE), # One-step fraction for pattern
    total_count = n(),
    .groups = 'drop'
  ) %>%
  # Filter to patterns with sufficient sample size (or selected representatives)
  filter(total_count > 100) %>% 
  arrange(onestep_fraction)

# Add global "All" average as reference
all_avg <- data.frame(
  mutation_pattern = "All",
  mechanism = "All",
  onestep_fraction = mean(mnv1000g_classified$is_onestep, na.rm = TRUE),
  total_count = nrow(mnv1000g_classified)
)

plot_final <- bind_rows(plot_pattern_df, all_avg) %>%
  mutate(mutation_pattern = factor(mutation_pattern, levels = unique(mutation_pattern)))

# --- 2. Scatter plot ---
ggplot(plot_final, aes(x = mutation_pattern, y = onestep_fraction, color = mechanism)) +
  geom_point(size = 3) +
  # Color scheme to match image_7f8ccf.png
  scale_color_manual(values = c(
    "Ti at CpG" = "#56B4E9", 
    "Pol-zeta error" = "#9970AB", 
    "Slippage at repeat" = "#8C510A", 
    "Tv combination" = "#D6604D"  )) +
  # Y-axis range
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "MNV pattern",
    y = "Fraction of one-step MNVs",
    color = "Mechanism"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )




library(dplyr)
library(tidyr)
library(ggplot2)




# Load TFBS prediction results (from file shown in image_824ae6.png)
# Ensure matrix_id (e.g., MA0073.1) matches JASPAR table
tfbs_res <- read.table("/home/caow/03mnv/analyse3/06TFBS/04res_mRNA/jaspar_pair.txt", header = FALSE, stringsAsFactors = FALSE)
# Based on the figure, assume column 1 is matrix_id and 4th from last is Effect (Gain/Loss)
colnames(tfbs_res)[c(1, 15)] <- c("matrix_id", "effect")

# Load JASPAR family metadata
jaspar_meta <- read.table("/home/caow/03mnv/analyse3/00addition/data/ultimate_metadata_table_CORE.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)


# Compare MNV/SNV rate in TF genes vs non-TF genes
tfgenelist <- jaspar_meta %>%
  transmute(name = toupper(name)) %>%
  separate_rows(name, sep = "::") %>%        # Split into multiple rows
  mutate(name = str_trim(name)) %>%          # Trim whitespace (safe)
  filter(name != "") %>%
  distinct(name)
colnames(tfgenelist)[1] <- 'gene'

snvnum1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000ggene_snv_counts.tsv',col_names = F) %>% filter(X2!='.')
mnvnum1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000ggene_mnv_counts.tsv',col_names = F) %>% filter(X2!='.')
allgenelist <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/hg38genes.bed',col_names = F) %>% filter(X5!='.')
colnames(snvnum1000g)[2] <- 'gene'
colnames(mnvnum1000g)[2] <- 'gene'
colnames(allgenelist)[5] <- 'gene'
allgenelist <- allgenelist %>% distinct(gene)
allgenelist %>% left_join(snvnum1000g,by = 'gene') %>% left_join(tfgenelist,by='gene')
mnvnum1000g %>% group_by(gene) %>% summarise(mnv_count=sum(X4))

final_result <- allgenelist %>%
  # 1. Flag TF genes
  # If gene is in tfgenelist$gene, mark TRUE; else FALSE
  mutate(is_tf = gene %in% tfgenelist$gene) %>%
  
  # 2. Join SNV counts from snvnum1000g (X4 column)
  # Select only gene and X4 to avoid extra columns
  left_join(snvnum1000g %>% select(gene, snv_count = X4), by = "gene") %>%
  
  # 3. If snvnum1000g missing (NA), set to 0
  mutate(snv_count = replace_na(snv_count, 0)) %>% 
  left_join(mnvnum1000g %>% group_by(gene) %>% summarise(mnv_count=sum(X4)), by = "gene") %>%
  
  # 3. If mnvnum1000g missing (NA), set to 0
  mutate(mnv_count = replace_na(mnv_count, 0))

plot_data <- final_result %>%
  filter(snv_count > 0) %>% # Drop genes without SNVs
  mutate(
    ratio = mnv_count / snv_count,
    Group = ifelse(is_tf, "TF Genes", "Non-TF Genes")
  )

# 2. Dynamically compute Y-axis max (auto-fit)
# Use max upper whisker from both groups as Y max
# Shows full boxes while trimming extreme outliers
# boxplot.stats(x)$stats[5] returns upper whisker (typically Q3 + 1.5*IQR)
ymax_tf <- boxplot.stats(plot_data$ratio[plot_data$Group == "TF Genes"])$stats[5]
ymax_non <- boxplot.stats(plot_data$ratio[plot_data$Group == "Non-TF Genes"])$stats[5]
y_limit_upper <- max(ymax_tf, ymax_non) * 1.2 # Multiply by 1.2 for P-value headroom

# 3. Plot
ggplot(plot_data, aes(x = Group, y = ratio, fill = Group)) +
  # --- Key change 1: hide outliers ---
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  
  # Add significance test (Wilcoxon)
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.5, 
                     label.y = y_limit_upper * 0.9, # Adjust P-value position
                     size = 5) +
  
  # Colors
  scale_fill_manual(values = c("TF Genes" = "#E41A1C", "Non-TF Genes" = "#377EB8")) +
  
  # --- Key change 2: adjust view range (Zoom In) ---
  # Use coord_cartesian to hide outliers without removing them
  coord_cartesian(ylim = c(0, y_limit_upper)) +
  
  labs(
    # title = "Comparison of MNV/SNV Ratio (Outliers Removed)",
    # subtitle = "Wilcoxon rank-sum test",
    y = "MNV/SNV Ratio",
    x = NULL
  ) +
  
  theme_bw(base_size = 14) + mythem+
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5)
  )

# Potential motif binding in ref
refback <- read_table('~/mnv2tfbs/data1118/hg38/all_motif_counts.txt',col_names = F)
colnames(refback) <- c('num','matrix_id')

classback <- jaspar_meta %>% select(matrix_id,class) %>% inner_join(refback,by= "matrix_id") %>% 
  group_by(class) %>% summarise(all=sum(num))

# 1. Join JASPAR class info
merged_data <- tfbs_res %>%
  inner_join(jaspar_meta %>% select(matrix_id, class, family), by = "matrix_id")

# 2. Count Gain and Loss per class
class_stats <- merged_data %>%
  group_by(class, effect) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = effect, values_from = count, values_fill = 0)

# 3. Compute rates
# Note: rates require a denominator; if using total predictions per class:
class_stats <- class_stats %>%
  mutate(Total = Gain + Loss,
         Gain_Rate = Gain / Total,
         Loss_Rate = Loss / Total) %>%
  arrange(desc(Total)) %>% inner_join(classback,by='class')

tf_data <- class_stats %>% filter(class!='Uncharacterized')
# --- 1. Data processing and normalization ---
plot_data <- tf_data %>%
  # Compute global baseline rates
  mutate(
    global_gain_rate = sum(Gain) / sum(all),
    global_loss_rate = sum(Loss) / sum(all)
  ) %>%
  rowwise() %>%
  mutate(
    # Compute raw rate per class
    raw_rate_gain = Gain / all,
    raw_rate_loss = Loss / all,
    
    # Enrichment (O/E): class rate / global average
    # >1 indicates enrichment; <1 indicates depletion
    Enrichment_Gain = raw_rate_gain / global_gain_rate,
    Enrichment_Loss = raw_rate_loss / global_loss_rate
  ) %>%
  ungroup() %>%
  # Pivot to long format for plotting
  select(class, Enrichment_Gain, Enrichment_Loss) %>%
  pivot_longer(cols = c(Enrichment_Gain, Enrichment_Loss), 
               names_to = "Type", values_to = "Enrichment") %>%
  mutate(
    Type = ifelse(Type == "Enrichment_Gain", "Gain (Creating)", "Loss (Disrupting)"),
    # Order classes by Gain enrichment
    class = factor(class, levels = tf_data$class[order(tf_data$Gain / tf_data$all)])
  )

# Recompute with Gain and Loss combined
# --- 1. Data processing and normalization ---
plot_data <- tf_data %>%
  # Compute global baseline rate
  mutate(
    global_rate = sum(Total) / sum(all),
  ) %>%
  rowwise() %>%
  mutate(
    # Compute raw rate per class
    raw_rate = Total / all,

    # Enrichment (O/E): class rate / global average
    # >1 indicates enrichment; <1 indicates depletion
    Enrichment_change = raw_rate / global_rate,
  ) %>%
  ungroup() %>%
  # Pivot to long format for plotting
  select(class, Enrichment_change) %>%
  pivot_longer(cols = c(Enrichment_change), 
               names_to = "Type", values_to = "Enrichment") %>%
  mutate(
    Type = 'Disrupted',
    # Order classes by enrichment
    class = factor(class, levels = tf_data$class[order(tf_data$Total / tf_data$all)])
  )
mythem <- theme(
  panel.grid = element_blank(),
  plot.margin = unit(rep(2, 4), 'lines'),
  panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
  text = element_text(size = 8, face = "bold"),
  axis.text.x = element_text(vjust = 1, size = 8, face = "bold", colour = 'black'),
  axis.title.x = element_text(vjust = 1, size = 8, face = "bold", color = 'black'),
  axis.text.y = element_text(hjust = 1, size = 8, face = "bold", color = 'black'),
  axis.title.y = element_text(vjust = 1, size = 8, face = "bold", color = 'black')
)

plot_data$log2FC <- log2(plot_data$Enrichment)

ggplot(plot_data, aes(x = log2FC, y = reorder(class, log2FC))) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_segment(aes(xend = 0, yend = class, color = log2FC > 0), size = 1.2) +
  geom_point(aes(color = log2FC > 0), size = 4) +
  
  scale_color_manual(values = c("FALSE" = "#4393C3", "TRUE" = "#D6604D"), 
                     labels = c("Depleted", "Enriched"), 
                     name = "Status") +
  
  labs(
    x = "Log2 (Enrichment Ratio)", # Update axis label if needed
    y = NULL
  ) +
  theme_bw()+mythem+theme(legend.position = "none")

# 1. Compute class order from log2FC (low to high or vice versa)
# Here we sort by log2FC ascending so it plots bottom to top
class_order <- plot_data$class[order(plot_data$log2FC)]

# 2. Set class as a factor with the same levels in both data frames
plot_data$class <- factor(plot_data$class, levels = class_order)

gc_data <- merged_data  %>% filter(class!='Uncharacterized') %>% 
  # --- Step A: compute GC content ---
  # Use stringr::str_count for G and C, divide by total length
  mutate(gc_content = str_count(V14, "[GC]") / nchar(V14))

# Note: gc_data is assumed to be your GC-content data frame
# Ensure gc_data$class matches plot_data$class
gc_data$class <- factor(gc_data$class, levels = class_order) 


# =======================================================
# Step 2: plot left lollipop chart (minor tweaks)
# =======================================================
p1 <- ggplot(plot_data, aes(x = log2FC, y = class)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(xend = 0, yend = class, color = log2FC > 0), size = 1.2) +
  geom_point(aes(color = log2FC > 0), size = 4) +
  scale_color_manual(values = c("FALSE" = "#4393C3", "TRUE" = "#D6604D"), 
                     labels = c("Depleted", "Enriched"), 
                     name = "Status") +
  labs(
    x = "Log2 (Enrichment Ratio)", 
    y = NULL
  ) +
  theme_bw() + 
  mythem + # Uncomment if you use a custom theme
  theme(legend.position = "none")


# =======================================================
# Step 3: plot right-side boxplot
# =======================================================
p2 <- ggplot(gc_data, aes(x = gc_content, y = class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "#E8F1F2") + # Box color is configurable
  labs(
    x = "GC Content", 
    y = NULL
  ) +
  theme_bw() + 
  mythem + # Keep theme consistent
  theme(
    # Key: hide Y-axis text/ticks/title because the left plot already shows them
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank() # Optional: remove horizontal grid
  )



# =======================================================
# Step 4: assemble plots
# =======================================================
# Use patchwork for layout
# widths controls relative widths, e.g. c(1.5, 1) makes left wider
# 1. Update left plot p1: remove right margin (r = 0)
# t=top, r=right, b=bottom, l=left
p1 <- p1 + theme(plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5, unit = "pt"))

# 2. Update right plot p2: remove left margin (l = 0)
p2 <- p2 + theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0, unit = "pt"))

# 3. Combine
combined_plot <- p1 + p2 + plot_layout(widths = c(1.5, 1.3))

# 1. Merge two data frames with aligned class values
# inner_join aligns by class column
correlation_data <- plot_data %>%
  select(class, log2FC) %>%
  inner_join(avg_gc_data, by = "class")
# Compute Spearman correlation
res <- cor.test(correlation_data$log2FC, 
                correlation_data$Mean_GC, 
                method = "spearman")

# Print full result
print(res)

# Extract key stats for downstream use
rho_value <- round(res$estimate, 3) # Correlation coefficient rho
p_value   <- formatC(res$p.value, format = "e", digits = 2) # P-value (scientific notation)

# One-line summary
paste0("Spearman Rho = ", rho_value, ", P-value = ", p_value)

p_corr <- ggplot(correlation_data, aes(x = log2FC, y = Mean_GC)) +
  geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
  geom_smooth(method = "lm", color = "red", se = FALSE, linetype = "dashed") + # Add fitted line
  
  # Auto-add correlation and P-value labels
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  
  labs(
    x = "Log2 (Enrichment Ratio)",
    y = "Average GC Content",
  ) +
  theme_bw()+mythem
avg_gc_data <- gc_data %>%
  group_by(class) %>%
  summarise(Mean_GC = mean(gc_content, na.rm = TRUE)) # Compute mean

# Key: ensure right-plot classes use the same factor order
avg_gc_data$class <- factor(avg_gc_data$class, levels = class_order)
p2 <- ggplot(avg_gc_data, aes(x = Mean_GC, y = class)) +
  geom_col(width = 0.6, fill = "#7F8C8D", alpha = 0.8) + # Neutral gray, slightly thinner bars
  
  # To show numeric labels beside bars, enable the next line:
  # geom_text(aes(label = round(Mean_GC, 2)), hjust = -0.2, size = 3) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) + # Tight to Y axis with right padding
  labs(x = "Average GC Content", y = NULL) +
  theme_bw() +
  theme(
    # Key: hide all Y-axis elements for seamless alignment
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(), # Remove horizontal grid for a cleaner look
    plot.margin = margin(l = 0) # Remove left margin to bring plots closer
  )


# =======================================================
# Step 4: combine
# =======================================================
# Adjust widths; typically Enrichment main, GC secondary (e.g., 2:1)
combined_plot <- p1 + p2 + plot_layout(widths = c(1.8, 1))
# --- 2. Lollipop chart ---
ggplot(plot_data, aes(x = Enrichment, y = class, color = Type)) +
  # Baseline x=1 (average)
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  
  # Connecting lines
  geom_line(aes(group = class), color = "grey80", size = 0.5) +
  
  # Points
  geom_point(size = 3.5) +
  
  # Colors (Gain red, Loss blue)
  scale_color_manual(values = c("Gain (Creating)" = "#D6604D", "Loss (Disrupting)" = "#4393C3")) +
  
  labs(
    title = "Relative Enrichment of MNV-induced TFBS Changes",
    subtitle = "Normalized by total motif count (O/E Ratio)",
    x = "Fold Enrichment (Observed Rate / Global Average)",
    y = NULL,
    color = "Event Type"
  ) +
  
  theme_bw() +mythem

# Prepare plot data (long format)
plot_data <- class_stats %>%
  select(class, Gain, Loss) %>%
  pivot_longer(cols = c(Gain, Loss), names_to = "Effect", values_to = "Count")



# Stacked bar chart showing impact across TF classes
ggplot(plot_data, aes(x = reorder(class, Count), y = Count, fill = Effect)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Horizontal layout for long class names
  scale_fill_manual(values = c("Gain" = "#D6604D", "Loss" = "#438fa9")) +
  labs(title = "Impact of MNVs on different TF Classes",
       x = "TF Class (JASPAR)",
       y = "Number of affected TFBS") +
  theme_bw()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          