# 各个数据集的区域分布情况
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
  
  
  # 读取SNV数量
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

# 读取gnomAD约束算出的各区域期望变异数
expected_var <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/expected_by_type.hg38.tsv')
# 1. 确保 expected_var 的顺序与 df 完全一致
# 使用 left_join 是最稳妥的方法，避免因 type 缺失或多余导致报错
df_merge <- df %>%
  left_join(expected_var, by = "type")

# 2. 定义需要参与计算的列（数据集列）
target_cols <- c("1000G", "GTEx", "UKB50w", "TCGA", "UKB20w")

# 3. 执行相除计算
# 我们创建一个新表来存储 ratio 结果
df_ratio <- df_merge

df_ratio[target_cols] <- sweep(df_merge[target_cols], 1, df_merge$expected_sum, "/")

# 4. 查看结果
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


# 绘制五条折线图
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

# --- 1. 主题设置 (保持不变) ---
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

# --- 2. 数据读取 ---
fun1 <- function(f) {
  df <- read.table(f, header = F, stringsAsFactors = F, sep = '\t')
  colnames(df) <- c('type', 'num')
  return(df)
}

# 读取 MNV 原始数量 (注意：这里先保留原始计数，不要覆盖)
df_mnv <- cbind(fun1('01num/1000G.txt'), fun1('01num/GTEx.txt'), fun1('01num/UKB50w.txt'), fun1('01num/TCGA.txt'), fun1('01num/UKB20w.txt'))
df_mnv <- df_mnv[c(1:2, 4, 6, 8, 10)]
colnames(df_mnv) <- c('type', '1000G', 'GTEx', 'UKB50w', 'TCGA', 'UKB20w')
df_mnv$id <- 1:nrow(df_mnv)

# 读取 SNV 原始数量
df_snv <- cbind(fun1('01num/snv_1000G.region.counts.tsv'), fun1('01num/snv_GTEx.region.counts.tsv'), fun1('01num/snv_UKB50w.region.counts.tsv'), fun1('01num/snv_TCGA.region.counts.tsv'), fun1('01num/snv_UKB20w.region.counts.tsv'))
df_snv <- df_snv[c(1:2, 4, 6, 8, 10)]
colnames(df_snv) <- c('type', '1000G', 'GTEx', 'UKB50w', 'TCGA', 'UKB20w')
df_snv$id <- 1:nrow(df_snv)

# --- 3. 计算 Ratio 数据框 ---
# 复制一份用于存储比率
df_ratio <- df_mnv
# 计算每个区域的 Ratio = MNV / SNV
df_ratio[, 2:6] <- df_mnv[, 2:6] / df_snv[, 2:6]



# --- 4. 定义绘图函数 (包含新的基准线计算逻辑) ---
plot_single_dataset <- function(data_ratio, counts_mnv, counts_snv, col_name, plot_title, line_color) {
  
  # 提取绘图用的比率数据
  plot_data <- data.frame(
    id = data_ratio$id,
    type = data_ratio$type,
    value = data_ratio[[col_name]]
  )
  
  # --- 关键修改：计算全局基准线 (Global Baseline) ---
  # 定义：该数据集中所有 MNV 之和 / 所有 SNV 之和
  # total_mnv <- sum(counts_mnv[[col_name]], na.rm = TRUE)
  # total_snv <- sum(counts_snv[[col_name]], na.rm = TRUE)
  # baseline_value <- 0.02269008830602263
  # 
  baseline_value <- mean(plot_data$value, na.rm = TRUE)
  
  print(paste("Baseline for", col_name, ":", baseline_value)) # 打印出来检查一下
  
  p <- ggplot(plot_data, aes(x = id, y = value)) +
    # 1. 柱状图 (带透明度)
    geom_bar(stat = "identity", fill = line_color, alpha = 0.3, width = 0.6) +
    
    # 2. 折线和点
    geom_line(color = line_color, lwd = 1) +
    geom_point(size = 2, color = line_color) +
    
    # 3. 基准虚线 (Global Rate)
    geom_hline(yintercept = baseline_value, linetype = "dashed", color = "grey30", size = 0.8) +
    
    # 4. 坐标轴设置
    labs(title = plot_title, y = 'MNV/SNV Ratio',x=NULL) +
    scale_x_continuous(breaks = c(1:nrow(data_ratio)), labels = gsub("_", " ", data_ratio$type)) +
    
    # 5. 主题
    theme_bw() + mythem
  
  return(p)
}

# --- 5. 执行绘图 ---

# 画 1000G
p_1000G <- plot_single_dataset(
  data_ratio = df_ratio,
  counts_mnv = df_mnv,   # 传入原始计数用于算总和
  counts_snv = df_snv,   # 传入原始计数用于算总和
  col_name = "1000G", 
  plot_title = "1000G", 
  line_color = "#d02f43"
)

# 画 GTEx
p_GTEx <- plot_single_dataset(
  data_ratio = df_ratio,
  counts_mnv = df_mnv, 
  counts_snv = df_snv,
  col_name = "GTEx", 
  plot_title = "GTEx", 
  line_color = "#e35e42"
)

# --- 6. 组合输出 ---
grid.arrange(p_1000G, p_GTEx, nrow = 1)

# 方法B: 分别保存为 PDF
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
## 读入数据
{
  mnv_diff<-read_tsv('02RNAFold/mnv/mnv_diff.res',col_names = F)
  snv_diff<-read_tsv('02RNAFold/snv/snv_diff.res',col_names = F)
  mnv_snv_diff<-read_tsv('02RNAFold/merge_res.txt',col_names = F)
}

# 2. 数据处理 (支持任意数量的 SNV 位点)
plot_data <- mnv_snv_diff %>%
  rowwise() %>% # 开启“按行处理”模式
  mutate(
    # 核心逻辑：
    # 1. str_split(X4, ","): 按逗号分割字符串
    # 2. unlist: 将分割后的列表解开成向量
    # 3. as.numeric: 转为数字
    # 4. sum: 求和
    sum_snv_effects = sum(as.numeric(unlist(str_split(X4, ",")))),
    
    # Y 轴：MNV effect
    mnv_effect = X3
  ) %>%
  ungroup() # 处理完后取消按行模式，提高后续运行效率

# 3. 定义你的主题 (保持不变)
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
# 4. 绘图
p <- ggplot(plot_data, aes(x = sum_snv_effects, y = mnv_effect)) +
  geom_point(color = "#377EB8", alpha = 0.7, size = 1.5) + 
  
  # 添加对角线 (y=x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#E41A1C", size = 0.8) +
  
  # 设置坐标轴标签
  # 注意：这里把标签改成更通用的写法，不再特指 SNV1+SNV2
  labs(
    x = expression(bold(Sum~of~individual~SNV~effects~(Sigma*Delta*MFE[SNV]))),
    y = expression(bold(MNV~effect~(Delta*MFE[MNV])))
  ) +
  
  theme_bw() + 
  mythem
p_with_cor <- p + 
  stat_cor(
    method = "pearson",      # 使用皮尔逊相关系数
    label.x.npc = 0.05,      # 标签水平位置 (0-1之间，0.05表示靠左)
    label.y.npc = 0.95,      # 标签垂直位置 (0-1之间，0.95表示靠上)
    # 或者你可以使用具体坐标，例如: label.x = -10, label.y = 10
    size = 4,                # 字体大小
    fontface = "bold"        # 字体加粗
  )
# 5. 显示并保存
# ggsave("MNV_vs_MultiSNVs_scatter.pdf", p, width = 5, height = 5)

threshold_value <- 1.0 

# 2. 计算偏差和比例
stats_result <- plot_data %>%
  mutate(
    # 计算实际值与预测值的差 (Residuals)
    deviation = mnv_effect - sum_snv_effects,
    
    # 取绝对值，判断是否超过阈值
    # abs() 是为了不管它是变得更稳定还是更不稳定，只要偏离够大就算
    is_non_additive = abs(deviation) > threshold_value
  )

# 3. 统计具体数值
total_count <- nrow(stats_result)
non_additive_count <- sum(stats_result$is_non_additive)
percentage <- round((non_additive_count / total_count) * 100, 2)

mnvcategories <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/mnv_categories_1000g_gtex')
mnv1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000G_mnv_annotation')
mnvgtex <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/GTEx_mnv_annotation')
mnv1000g_anno <- mnv1000g %>% inner_join(mnvcategories,by='MNVid')
mnvgtex_anno <- mnvgtex %>% inner_join(mnvcategories,by='MNVid')

custom_colors <- c(
  "snv_event" = "#A83232",   # 深砖红色
  "one_step"  = "#E8A03D",   # 橙色
  "repeat"    = "#2C4B8E"    # 深蓝色
)

# B. 定义标签文字 (图2中的文字)
custom_labels <- c(
  "snv_event" = "SNV event",
  "one_step"  = "One-step",
  "repeat"    = "Slippage at Repeat"
)

my_order <- c('UTR5','exon','splice','intron','UTR3','lncRNA','miRNA','circRNA','piRNA','ATAC','conserved_element','enhancer','miRBS','TFBS','intergenic')
# 2. 数据清洗与转换 (核心步骤)
plot_data <- mnv1000g_anno %>%
  # 第一步：将含有分号的 annotation 拆分成多行
  # 例如 "lncRNA;piRNA" 会变成两行，分别对应 lncRNA 和 piRNA
  separate_rows(annotation, sep = ";") %>%
  
  # 第二步：将三列分类转换为长格式
  # 我们只需要值为 1 的行，因为那代表该 MNV 属于该类别
  pivot_longer(
    cols = c(snv_event, one_step, `repeat`), 
    names_to = "category", 
    values_to = "is_present"
  ) %>%
  
  # 过滤掉值为0的行（即不属于该类别的行）
  filter(is_present == 1)

plot_data <- plot_data %>%
  mutate(annotation = factor(annotation, levels = my_order))

my_order <- c('UTR5','exon','splice','intron','UTR3','lncRNA','miRNA','circRNA','piRNA','ATAC','conserved_element','enhancer','miRBS','TFBS','intergenic','TSS')
# 2. 数据清洗与转换 (核心步骤)
plot_data <- mnv1000g_anno %>%
  # 第一步：将含有分号的 annotation 拆分成多行
  # 例如 "lncRNA;piRNA" 会变成两行，分别对应 lncRNA 和 piRNA
  separate_rows(annotation, sep = ";") %>%
  
  # 第二步：将三列分类转换为长格式
  # 我们只需要值为 1 的行，因为那代表该 MNV 属于该类别
  pivot_longer(
    cols = c(snv_event, one_step, `repeat`), 
    names_to = "category", 
    values_to = "is_present"
  ) %>%
  
  # 过滤掉值为0的行（即不属于该类别的行）
  filter(is_present == 1)

# 比较了不同类别MNV在各个基因组功能区域的比例
plot_data <- plot_data %>%
  mutate(annotation = factor(annotation, levels = my_order))

mnvtypenum <- plot_data %>% group_by(annotation,category) %>% summarise(n=n())
colnames(mnvtypenum) <- c('type','category','1000G')

plot_data1_gtex <- mnvgtex_anno %>% distinct() %>% 
  # 第一步：将含有分号的 annotation 拆分成多行
  # 例如 "lncRNA;piRNA" 会变成两行，分别对应 lncRNA 和 piRNA
  separate_rows(annotation, sep = ";") %>%
  
  # 第二步：将三列分类转换为长格式
  # 我们只需要值为 1 的行，因为那代表该 MNV 属于该类别
  pivot_longer(
    cols = c(snv_event, one_step, `repeat`), 
    names_to = "category", 
    values_to = "is_present"
  ) %>%
  
  # 过滤掉值为0的行（即不属于该类别的行）
  filter(is_present == 1)

plot_data1_gtex <- plot_data1_gtex %>%
  mutate(annotation = factor(annotation, levels = my_order))

mnvtypenum_gtex <- plot_data1_gtex %>% group_by(annotation,category) %>% summarise(n=n())
colnames(mnvtypenum_gtex) <- c('type','category','GTEx')

mnvtypenum <-  mnvtypenum %>% left_join(mnvtypenum_gtex,by=c('type','category'))
mnvtypenum1 <-  mnvtypenum %>% filter(category=='one_step') %>% select(-2)
mnvtypenum2 <-  mnvtypenum %>% filter(category=='repeat') %>% select(-2)
mnvtypenum3 <-  mnvtypenum %>% filter(category=='snv_event') %>% select(-2)


# 读取 SNV 原始数量
df_snv <- cbind(fun1('01num/snv_1000G.region.counts.tsv'), fun1('01num/snv_GTEx.region.counts.tsv'), fun1('01num/snv_UKB50w.region.counts.tsv'), fun1('01num/snv_TCGA.region.counts.tsv'), fun1('01num/snv_UKB20w.region.counts.tsv'))
df_snv <- df_snv[c(1:2, 4, 6, 8, 10)]
colnames(df_snv) <- c('type', '1000G', 'GTEx', 'UKB50w', 'TCGA', 'UKB20w')
df_snv$id <- 1:nrow(df_snv)

# --- 3. 计算 Ratio 数据框 ---
# 复制一份用于存储比率
# df_mnv <- df_mnv[-17,]

df_ratio1 <- df_mnv
df_ratio2 <- df_mnv
df_ratio3 <- df_mnv

# df_snv <- df_snv[-16,]

# 计算每个区域的 Ratio = MNV / SNV
df_ratio1[, 2:3] <- mnvtypenum1[-17, 2:3] / df_snv[, 2:3]
df_ratio2[, 2:3] <- mnvtypenum2[-17, 2:3] / df_snv[, 2:3]
df_ratio3[, 2:3] <- mnvtypenum3[-17, 2:3] / df_snv[, 2:3]

plot_single_dataset1 <- function(data_ratio, counts_mnv, counts_snv, col_name, plot_title, line_color, baseline) {
  
  # 提取绘图用的比率数据
  plot_data <- data.frame(
    id = data_ratio$id,
    type = data_ratio$type,
    value = data_ratio[[col_name]]
  )
  
  # --- 关键修改：计算全局基准线 (Global Baseline) ---
  # 定义：该数据集中所有 MNV 之和 / 所有 SNV 之和
  # total_mnv <- sum(counts_mnv[[col_name]], na.rm = TRUE)
  # total_snv <- sum(counts_snv[[col_name]], na.rm = TRUE)
  # baseline_value <- 0.02269008830602263
  # 
  baseline_value <- baseline
  
  print(paste("Baseline for", col_name, ":", baseline_value)) # 打印出来检查一下
  
  p <- ggplot(plot_data, aes(x = id, y = value)) +
    # 1. 柱状图 (带透明度)
    geom_bar(stat = "identity", fill = line_color, alpha = 0.3, width = 0.6) +
    
    # 2. 折线和点
    geom_line(color = line_color, lwd = 1) +
    geom_point(size = 2, color = line_color) +
    
    # 3. 基准虚线 (Global Rate)
    geom_hline(yintercept = baseline_value, linetype = "dashed", color = "grey30", size = 0.8) +
    
    # 4. 坐标轴设置
    labs(title = plot_title, x = 'Genomic Region', y = 'MNV-to-SNV Ratio') +
    scale_x_continuous(breaks = c(1:nrow(data_ratio)), labels = data_ratio$type) +
    
    # 5. 主题
    theme_bw() + mythem
  
  return(p)
}

# 画 1000G
p_1000G1 <- plot_single_dataset1(
  data_ratio = df_ratio1,
  counts_mnv = mnvtypenum1,   # 传入原始计数用于算总和
  counts_snv = df_snv,   # 传入原始计数用于算总和
  col_name = "1000G", 
  plot_title = "1000G Dataset", 
  line_color = "#d02f43",baseline = 0.041212660155921
)

# 画 GTEx
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
  counts_mnv = mnvtypenum1,   # 传入原始计数用于算总和
  counts_snv = df_snv,   # 传入原始计数用于算总和
  col_name = "1000G", 
  plot_title = "1000G Dataset", 
  line_color = "#d02f43",baseline = 0.041212660155921
)

# 画 GTEx
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
  counts_mnv = mnvtypenum1,   # 传入原始计数用于算总和
  counts_snv = df_snv,   # 传入原始计数用于算总和
  col_name = "1000G", 
  plot_title = "1000G Dataset", 
  line_color = "#d02f43",baseline = 0.041212660155921
)

# 画 GTEx
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

# 定义每一行的标题样式
title_style <- function(label) {
  textGrob(label, gp = gpar(fontsize = 12, fontface = "bold"))
}

# 第一步：将每一行与其对应的机制标签组合
# Row 1: SNV event
row1 <- arrangeGrob(p_1000G1, p_GTEx1, ncol = 2, 
                    top = title_style("pol-zeta-mediated one-step"))

# Row 2: pol-zeta-mediated one-step
row2 <- arrangeGrob(p_1000G2, p_GTEx2, ncol = 2, 
                    top = title_style("repeat slippage"))

# Row 3: repeat slippage
row3 <- arrangeGrob(p_1000G3, p_GTEx3, ncol = 2, 
                    top = title_style("SNV event"))

# 第二步：将三行垂直排列
grid.arrange(row1, row2, row3, nrow = 3)


prepare_plot_data <- function(df1, df2, df3, dataset_col) {
  # 提取 Pol-zeta
  d1 <- df1[, c("id", "type", dataset_col)]
  colnames(d1)[3] <- "value"
  d1$Type <- "Pol-zeta-mediated one-step"
  
  # 提取 Repeat slippage
  d2 <- df2[, c("id", "type", dataset_col)]
  colnames(d2)[3] <- "value"
  d2$Type <- "Repeat slippage"
  
  # 提取 SNV event
  d3 <- df3[, c("id", "type", dataset_col)]
  colnames(d3)[3] <- "value"
  d3$Type <- "SNV event"
  
  # 合并数据
  plot_data <- rbind(d1, d2, d3)
  
  # 确保 X 轴顺序：将 type 转为因子，顺序按 id 排序
  # 假设 id 是 1, 2, 3... 对应原本的行顺序
  order_levels <- unique(d1$type[order(d1$id)])
  plot_data$type <- factor(plot_data$type, levels = order_levels)
  
  # 确保 Mechanism 的图例顺序
  plot_data$Type <- factor(plot_data$Type, 
                                levels = c("SNV event", "Pol-zeta-mediated one-step", "Repeat slippage"))
  
  return(plot_data)
}

# =======================================================
# 2. 定义绘图函数 (多机制合一)
# =======================================================
plot_lines_final <- function(plot_data, plot_title, baseline_value, show_legend = TRUE, show_x_title = TRUE) {
  
  print(paste("Drawing:", plot_title, "| Baseline:", baseline_value)) 
  
  # 基础绘图对象
  p <- ggplot(plot_data, aes(x = type, y = value, color = Type, group = Type)) +
    
    # 1. 基准线
    # geom_hline(yintercept = baseline_value, linetype = "dashed", color = "grey40", size = 0.8) +
    
    # 2. 折线和点
    geom_line(lwd=1, alpha = 0.9) + 
    geom_point(size = 2) +

    # 3. 颜色设置
    scale_color_manual(values = c("SNV event" = "#A83232", 
                                  "Pol-zeta-mediated one-step" = "#E8A03D", 
                                  "Repeat slippage" = "#2C4B8E")) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    # 4. 标题与轴标签控制
    labs(title = plot_title, 
         # 根据参数决定是否显示 X 轴标题
         x = if(show_x_title) 'Genomic Region' else NULL, 
         y = 'MNV/SNV Ratio') +
    

    # 5. 主题设置
    theme_bw() + mythem+theme(legend.position = if(show_legend) "top" else "none")
    # theme(
    #   # --- 关键修改：标题居中 ---
    #   plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    #   
    #   # 坐标轴文字
    #   axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    #   axis.title.y = element_text(size = 12, face = "bold"),
    #   
    #   # --- 关键修改：图例位置控制 ---
    #   # 如果 show_legend 为 TRUE 放顶部，否则隐藏
    #   legend.position = if(show_legend) "top" else "none",
    #   legend.title = element_blank(),
    #   legend.text = element_text(size = 11, face = "bold"),
    #   
    #   panel.grid.minor = element_blank()
    #   # + mythem (如果需要)
    # )
  
  return(p)
}

# =======================================================
# 3. 执行绘图
# =======================================================

# --- 准备 1000G 数据 ---
# 注意：这里假设 df_ratio1, df_ratio2, df_ratio3 是你已经算好的三个dataframe
data_1000G <- prepare_plot_data(df_ratio1, df_ratio2, df_ratio3, "1000G")

# --- 准备 GTEx 数据 ---
data_GTEx <- prepare_plot_data(df_ratio1, df_ratio2, df_ratio3, "GTEx")

library(ggpubr) # install.packages("ggpubr")

# 1. 确保两个图的代码里 show_legend = TRUE 
# (ggpubr 会自动提取其中一个，并隐藏其他的)
p_1000G_with_leg <- plot_lines_final(data_1000G, "1000G", 0.0412, show_legend = TRUE, show_x_title = FALSE)
p_GTEx_with_leg  <- plot_lines_final(data_GTEx,  "GTEx",  0.0278, show_legend = TRUE, show_x_title = FALSE)

# 2. 一键排列 + 共享图例
p2 <- ggarrange(
  p_1000G_with_leg, 
  p_GTEx_with_leg, 
  ncol = 2, 
  nrow = 1,
  common.legend = TRUE, # 开启共享图例功能
  legend = "bottom"     # 图例放在最下面
)
p1 <- grid.arrange(p_1000G, p_GTEx, nrow = 1)
grid.arrange(p1, p2, ncol = 1)

library(cowplot)
library(ggpubr)

# --- 第一步：构建下方带有共享图例的图 (Row 2) ---
# 使用 ggarrange 提取共享图例
p2 <- ggarrange(
  p_1000G_with_leg, 
  p_GTEx_with_leg, 
  ncol = 2, 
  nrow = 1,
  common.legend = TRUE, 
  legend = "bottom"
)

# --- 第二步：构建上方的图 (Row 1) ---
# 这里也建议用 plot_grid 替代 grid.arrange，因为 plot_grid 返回的是 ggplot 对象，更稳定
p1 <- plot_grid(
  p_1000G, 
  p_GTEx, 
  ncol = 2, 
  align = "h" # 水平对齐
)

# --- 第三步：最终组合与间距调整 ---
final_plot <- plot_grid(
  p1, p2, 
  ncol = 1, 
  
  # 1. 【关键】对齐上下两行的坐标轴
  align = "v", 
  axis = "lr",
  
  # 2. 【关键】调整上下高度比例
  # 下面的图带图例和X轴标题，所以给它多一点空间 (例如 1:1.2)
  rel_heights = c(1, 1.1),
  
  # 3. 【关键】增加子图之间的物理间距 (可选)
  # scale < 1 会让图片稍微缩小，从而在图片周围产生留白(间距)
  scale = c(1.07, 1.07) 
)

final_plot

# --- 画图 ---
# 这里的 baseline 如果你想加，可以传参，或者在函数里写死
# 1000G 图 (Baseline: 0.0412...)
p_top <- plot_lines_final(
  plot_data = data_1000G, 
  plot_title = "1000G Dataset", 
  baseline_value = 0.041212660155921,
  show_legend = TRUE,     # 显示图例
  show_x_title = TRUE    # 隐藏 X 轴标题
)

# --- 下面的图 (GTEx) ---
# 去掉图例 (show_legend = FALSE) 因为上面已经有了
# 保留 X 轴标题 (show_x_title = TRUE)
p_bottom <- plot_lines_final(
  plot_data = data_GTEx, 
  plot_title = "GTEx Dataset", 
  baseline_value = 0.027862315966804,
  show_legend = FALSE,    # 隐藏图例
  show_x_title = TRUE     # 显示 X 轴标题
)

grid.arrange(p_top, p_bottom, nrow = 1)


library(cowplot) # 强烈建议使用这个包来拼图

plot_grid(p_top, p_bottom, 
          ncol = 1, 
          align = "h", # 垂直对齐坐标轴
          rel_heights = c(1.2, 1.1))


# 3. 绘图 (百分比堆叠图)
ggplot(plot_data, aes(x = annotation, fill = category)) +
  # position = "fill" 是关键，它会自动计算百分比并堆叠
  geom_bar(position = "stack", width = 0.7) + 
  
  # 添加百分比标签 (可选，为了更清晰)
  geom_text(stat = "count", 
            aes(label = scales::percent(after_stat(count) / tapply(after_stat(count), after_stat(x), sum)[after_stat(x)], accuracy = 1)),
            position = position_fill(vjust = 0.5), 
            size = 3, color = "black") +
  
  # 调整坐标轴和配色
  scale_y_continuous(labels = scales::percent) + # Y轴显示为百分比
  scale_fill_manual(
    values = custom_colors, # 应用我们定义的颜色向量
    labels = custom_labels, # 应用我们定义的标签向量
    name = NULL             # 去掉图例标题，让图例更干净 (或者写 name = "Variant Type")
  ) +
  labs(
    title = "Proportion of MNV Categories by Annotation Region",
    x = "Genomic Region (Annotation)",
    y = "Proportion",
    fill = "MNV Type"
  ) +
  theme_bw() + mythem+
  # 如果横坐标标签太长，可以让它们倾斜
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# 如果你想看具体的计数而不是比例，只需将 position = "fill" 改为 position = "stack"


mnv1000g_pattern <- read_tsv('~/03mnv/analyse3/00addition/data/mnv1000g_pattern',col_names = F)
colnames(mnv1000g_pattern) <- c('MNVid','ref','alt','dis')

# 创建 MNV 分类数据框
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

# 添加 mutation 列 (格式如 AA->CC)
mnv_data$mutation <- paste(mnv_data$refs, "->", mnv_data$alts, sep="")

# 调整列顺序，使其与图片更一致
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

# (可选) 如果你也要 r / r1 / r_all 的结构，这里可以建 list：
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
    # 比较 i 与 rc 在位置 1,3,5,7 的碱基（等价于 Python 的 range(0,8,2) 对字符串取 i[k]）
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
    # 理论上一定会 decided；若完全相等则其实 i==rc，前面已排除
  } else {
    mapdict[i] <- i
  }
}

# mapdict 就是最终结果
# 示例：查看前几个映射
head(mapdict)

mnv1000g_2joint_pattern <- mnv1000g %>% inner_join(mnv1000g_pattern,by='MNVid') %>% filter(nchar(ref)==3&dis==1)

# --- 第一步：构建全量对称查找表 ---
# 您的 mapdict 是 {smaller: larger} 或 {self: self}
# 我们需要一个表能把 larger 映射回 smaller，且 smaller 映射到自身
smaller_keys <- names(mapdict)
larger_values <- unname(mapdict)

# 构建一个全量查找向量：包含所有可能的原始模式到代表性模式的映射
full_lookup_vec <- setNames(smaller_keys, smaller_keys) # smaller -> smaller
full_lookup_vec <- c(full_lookup_vec, setNames(smaller_keys, larger_values)) # larger -> smaller
# 去重（防止自互补模式重复添加）
full_lookup_vec <- full_lookup_vec[!duplicated(names(full_lookup_vec))]

# 给mnv1000g添加prompoter和TSS
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

# --- 第二步：处理您的 mnv1000g 数据 ---
mnv1000g_classified <- mnv1000g %>% 
  inner_join(mnv1000g_pattern, by = 'MNVid') %>% 
  # 筛选 2-joint 且相邻的 MNV
  filter(nchar(ref) == 3 & dis == 1) %>% 
  mutate(
    # 1. 组合成 mapdict 的原始格式，例如 "C,G:T,A"
    raw_key = paste0(ref, ":", alt)
  ) %>% 
  mutate(
    # 2. 通过查找表直接获取代表性模式 (Representative Pattern)
    # 这步非常快，直接向量化替换
    rep_key = full_lookup_vec[raw_key]
  ) %>% 
  mutate(
    # 3. 将 "A,A:C,C" 转换为 "AA->CC" 格式以匹配 mnv_data
    match_pattern = str_replace_all(rep_key, ",", "") %>% str_replace(":", "->")
  ) %>% 
  # 4. 关联 mnv_data 获取 class
  left_join(mnv_data, by = c("match_pattern" = "mutation"))

mnvgtex_pattern <- read_tsv('~/03mnv/analyse3/00addition/data/mnvGTEx_pattern',col_names = F)
colnames(mnvgtex_pattern) <- c('MNVid','ref','alt','dis')
mnvgtex_2joint_pattern <- mnvgtex %>% inner_join(mnvgtex_pattern,by='MNVid') %>% filter(nchar(ref)==3&dis==1)

mnvgtex_classified <- mnvgtex %>% 
  inner_join(mnvgtex_pattern, by = 'MNVid') %>% 
  # 筛选 2-joint 且相邻的 MNV
  filter(nchar(ref) == 3 & dis == 1) %>% 
  mutate(
    # 1. 组合成 mapdict 的原始格式，例如 "C,G:T,A"
    raw_key = paste0(ref, ":", alt)
  ) %>% 
  mutate(
    # 2. 通过查找表直接获取代表性模式 (Representative Pattern)
    # 这步非常快，直接向量化替换
    rep_key = full_lookup_vec[raw_key]
  ) %>% 
  mutate(
    # 3. 将 "A,A:C,C" 转换为 "AA->CC" 格式以匹配 mnv_data
    match_pattern = str_replace_all(rep_key, ",", "") %>% str_replace(":", "->")
  ) %>% 
  # 4. 关联 mnv_data 获取 class
  left_join(mnv_data, by = c("match_pattern" = "mutation"))

# --- 第三步：查看统计结果 ---
# 统计各分类的数量
class_summary <- mnv1000g_classified %>%
  group_by(class) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

print(class_summary)

# 1. 定义突变模式列表 (根据你的要求)
cpg_list     <- c("GA->AG", "CC->TT", "AC->GT", "CA->TG")
polzeta_list <- c("GA->TT", "GC->AA")
rep_list     <- c("AC->CA", "CA->AC", "AT->TA", "TA->AT", "AA->TT", "AA->CC")
tv_list      <- c("TA->CG", "CG->AT", "AT->CG", "CG->GC", "GC->CG", "CG->AC")

# 2. 定义区域顺序
region_order <- c("enhancer", "Promoter", "TSS", "UTR5", "exon", "intron", "UTR3")

# 3. 数据处理
plot_df <- mnv1000g_classified %>%
  # 处理多重注释：按分号展开
  separate_rows(annotation, sep = ";") %>%
  
  # --- 核心修改：基于具体的 Pattern 映射主要突变机制 ---
  # 注意：这里使用的是最终匹配模式 (final_pattern 或 match_pattern)，格式为 "AA->CC"
  mutate(mechanism = case_when(
    match_pattern %in% cpg_list     ~ "Ti at CpG",
    match_pattern %in% polzeta_list ~ "Pol-zeta error",
    match_pattern %in% rep_list     ~ "Slippage at repeat",
    match_pattern %in% tv_list      ~ "Tv combination"
  )) %>%
  
  # 过滤目标区域
  filter(annotation %in% region_order) %>%
  
  # 按区域和机制统计比例
  group_by(annotation, mechanism) %>%
  summarise(count = n(), .groups = 'drop') %>%
  
  group_by(annotation) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup() %>%
  
  # 设置因子顺序
  mutate(annotation = factor(annotation, levels = region_order)) %>% filter(!is.na(mechanism))

# 4. 绘图
ggplot(plot_df, aes(x = annotation, y = fraction, color = mechanism, group = mechanism)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +
  
  # 配色方案 (保持与图片一致)
  scale_color_manual(values = c(
    "Ti at CpG"          = "#56B4E9", # 浅蓝
    "Pol-zeta error"     = "#9970AB", # 紫色
    "Slippage at repeat" = "#8C510A", # 棕色
    "Tv combination"     = "#D6604D" # 红色
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
  # 处理多重注释：按分号展开
  separate_rows(annotation, sep = ";") %>%
  
  # --- 核心修改：基于具体的 Pattern 映射主要突变机制 ---
  # 注意：这里使用的是最终匹配模式 (final_pattern 或 match_pattern)，格式为 "AA->CC"
  mutate(mechanism = case_when(
    match_pattern %in% cpg_list     ~ "Ti at CpG",
    match_pattern %in% polzeta_list ~ "Pol-zeta error",
    match_pattern %in% rep_list     ~ "Slippage at repeat",
    match_pattern %in% tv_list      ~ "Tv combination"
  )) %>%
  
  # 过滤目标区域
  filter(annotation %in% region_order) %>%
  
  # 按区域和机制统计比例
  group_by(annotation, mechanism) %>%
  summarise(count = n(), .groups = 'drop') %>%
  
  group_by(annotation) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup() %>%
  
  # 设置因子顺序
  mutate(annotation = factor(annotation, levels = region_order)) %>% filter(!is.na(mechanism))

# 4. 绘图
ggplot(plot_df1, aes(x = annotation, y = fraction, color = mechanism, group = mechanism)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +
  
  # 配色方案 (保持与图片一致)
  scale_color_manual(values = c(
    "Ti at CpG"          = "#56B4E9", # 浅蓝
    "Pol-zeta error"     = "#9970AB", # 紫色
    "Slippage at repeat" = "#8C510A", # 棕色
    "Tv combination"     = "#D6604D" # 红色
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
  # 筛选 2-joint 且相邻的 MNV
  filter(nchar(ref) == 3 & dis == 1) %>% 
  mutate(
    # 1. 组合成 mapdict 的原始格式，例如 "C,G:T,A"
    raw_key = paste0(ref, ":", alt)
  ) %>% 
  mutate(
    # 2. 通过查找表直接获取代表性模式 (Representative Pattern)
    # 这步非常快，直接向量化替换
    rep_key = full_lookup_vec[raw_key]
  ) %>% 
  mutate(
    # 3. 将 "A,A:C,C" 转换为 "AA->CC" 格式以匹配 mnv_data
    match_pattern = str_replace_all(rep_key, ",", "") %>% str_replace(":", "->")
  ) %>% 
  # 4. 关联 mnv_data 获取 class
  left_join(mnv_data, by = c("match_pattern" = "mutation"))
plot_df2 <- mnvmerge_classified %>%
  # 处理多重注释：按分号展开
  separate_rows(annotation, sep = ";") %>% distinct() %>% 
  
  # --- 核心修改：基于具体的 Pattern 映射主要突变机制 ---
  # 注意：这里使用的是最终匹配模式 (final_pattern 或 match_pattern)，格式为 "AA->CC"
  mutate(mechanism = case_when(
    match_pattern %in% cpg_list     ~ "Ti at CpG",
    match_pattern %in% polzeta_list ~ "Pol-zeta error",
    match_pattern %in% rep_list     ~ "Slippage at repeat",
    match_pattern %in% tv_list      ~ "Tv combination"
  )) %>%
  
  # 过滤目标区域
  filter(annotation %in% region_order) %>%
  
  # 按区域和机制统计比例
  group_by(annotation, mechanism) %>%
  summarise(count = n(), .groups = 'drop') %>%
  
  group_by(annotation) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup() %>%
  
  # 设置因子顺序
  mutate(annotation = factor(annotation, levels = region_order)) %>% filter(!is.na(mechanism))

# 4. 绘图
ggplot(plot_df2, aes(x = annotation, y = fraction, color = mechanism, group = mechanism)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +
  
  # 配色方案 (保持与图片一致)
  scale_color_manual(values = c(
    "Ti at CpG"          = "#56B4E9", # 浅蓝
    "Pol-zeta error"     = "#9970AB", # 紫色
    "Slippage at repeat" = "#8C510A", # 棕色
    "Tv combination"     = "#D6604D" # 红色
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



# --- 1. 准备绘图数据 ---
# 假设 mnv1000g_classified 中已有 'mechanism' 和判定 one-step 的逻辑
# 我们需要计算每个 pattern 的 one-step 比例

plot_pattern_df <- mnv1000g_classified %>%
  # 判定 one-step (此处需根据您原始数据中的 AC 列计算，示例逻辑如下)
  # mutate(is_onestep = ifelse(AC_mnv / AC_snv_min > 0.9 & (AC_snv_max - AC_snv_min)/AC_snv_min < 0.1, TRUE, FALSE)) %>%
  
  group_by(mutation_pattern, mechanism) %>%
  summarise(
    onestep_fraction = mean(is_onestep, na.rm = TRUE), # 计算该模式下 one-step 的比例
    total_count = n(),
    .groups = 'drop'
  ) %>%
  # 筛选样本量足够的模式进行展示，或选择图片中指定的代表性模式
  filter(total_count > 100) %>% 
  arrange(onestep_fraction)

# 添加一个 "All" 的全局平均值作为对照
all_avg <- data.frame(
  mutation_pattern = "All",
  mechanism = "All",
  onestep_fraction = mean(mnv1000g_classified$is_onestep, na.rm = TRUE),
  total_count = nrow(mnv1000g_classified)
)

plot_final <- bind_rows(plot_pattern_df, all_avg) %>%
  mutate(mutation_pattern = factor(mutation_pattern, levels = unique(mutation_pattern)))

# --- 2. 绘制散点图 ---
ggplot(plot_final, aes(x = mutation_pattern, y = onestep_fraction, color = mechanism)) +
  geom_point(size = 3) +
  # 设定颜色方案与 image_7f8ccf.png 一致
  scale_color_manual(values = c(
    "Ti at CpG" = "#56B4E9", 
    "Pol-zeta error" = "#9970AB", 
    "Slippage at repeat" = "#8C510A", 
    "Tv combination" = "#D6604D"  )) +
  # 设置 Y 轴范围
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




# 读取你的 TFBS 预测结果 (来自 image_824ae6.png 所示文件)
# 注意：确保 matrix_id 列（如 MA0073.1）与 JASPAR 表一致
tfbs_res <- read.table("/home/caow/03mnv/analyse3/06TFBS/04res_mRNA/jaspar_pair.txt", header = FALSE, stringsAsFactors = FALSE)
# 根据图片内容，假设第一列是 matrix_id，倒数第四列是 Effect (Gain/Loss)
colnames(tfbs_res)[c(1, 15)] <- c("matrix_id", "effect")

# 读取 JASPAR 家族信息表
jaspar_meta <- read.table("/home/caow/03mnv/analyse3/00addition/data/ultimate_metadata_table_CORE.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)


# 统计TF基因上MNV的发生率：MNV/SNV，与非TF基因比较
tfgenelist <- jaspar_meta %>%
  transmute(name = toupper(name)) %>%
  separate_rows(name, sep = "::") %>%        # 拆成多行
  mutate(name = str_trim(name)) %>%          # 去空格（保险）
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
  # 1. 统计是否为 TF 基因
  # 如果 gene 在 tfgenelist$gene 中，标记为 TRUE，否则为 FALSE
  mutate(is_tf = gene %in% tfgenelist$gene) %>%
  
  # 2. 合并 snvnum1000g 中的 SNV 数量 (X4列)
  # 我们先从 snvnum1000g 中只选出 gene 和 X4 两列，避免合并多余的列
  left_join(snvnum1000g %>% select(gene, snv_count = X4), by = "gene") %>%
  
  # 3. 如果在 snvnum1000g 找不到 (值为 NA)，则计为 0
  mutate(snv_count = replace_na(snv_count, 0)) %>% 
  left_join(mnvnum1000g %>% group_by(gene) %>% summarise(mnv_count=sum(X4)), by = "gene") %>%
  
  # 3. 如果在 snvnum1000g 找不到 (值为 NA)，则计为 0
  mutate(mnv_count = replace_na(mnv_count, 0))

plot_data <- final_result %>%
  filter(snv_count > 0) %>% # 过滤掉无 SNV 的基因
  mutate(
    ratio = mnv_count / snv_count,
    Group = ifelse(is_tf, "TF Genes", "Non-TF Genes")
  )

# 2. 动态计算 Y 轴上限 (为了自动适配数据)
# 我们取两组数据中 "箱线图上边缘 (Upper Whisker)" 的最大值作为 Y 轴上限
# 这样既能展示完整箱体，又截断了极端的离群值
# boxplot.stats(x)$stats[5] 返回的是上须的位置 (通常是 Q3 + 1.5*IQR)
ymax_tf <- boxplot.stats(plot_data$ratio[plot_data$Group == "TF Genes"])$stats[5]
ymax_non <- boxplot.stats(plot_data$ratio[plot_data$Group == "Non-TF Genes"])$stats[5]
y_limit_upper <- max(ymax_tf, ymax_non) * 1.2 # 乘以 1.2 留出一点顶部空间给 P 值

# 3. 绘图
ggplot(plot_data, aes(x = Group, y = ratio, fill = Group)) +
  # --- 核心修改 1: 隐藏离群点 ---
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  
  # 添加显著性检验 (Wilcoxon test)
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.5, 
                     label.y = y_limit_upper * 0.9, # 调整 P 值位置
                     size = 5) +
  
  # 配色
  scale_fill_manual(values = c("TF Genes" = "#E41A1C", "Non-TF Genes" = "#377EB8")) +
  
  # --- 核心修改 2: 调整视野范围 (Zoom In) ---
  # 使用 coord_cartesian 确保只是"看不见"离群值，而不是"删除"它们
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

# ref中潜在的motif结合
refback <- read_table('~/mnv2tfbs/data1118/hg38/all_motif_counts.txt',col_names = F)
colnames(refback) <- c('num','matrix_id')

classback <- jaspar_meta %>% select(matrix_id,class) %>% inner_join(refback,by= "matrix_id") %>% 
  group_by(class) %>% summarise(all=sum(num))

# 1. 关联 JASPAR 的 class 信息
merged_data <- tfbs_res %>%
  inner_join(jaspar_meta %>% select(matrix_id, class, family), by = "matrix_id")

# 2. 统计每个家族 (class) 中 Gain 和 Loss 的数量
class_stats <- merged_data %>%
  group_by(class, effect) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = effect, values_from = count, values_fill = 0)

# 3. 计算发生率 (Rate)
# 注意：发生率通常需要分母。如果你的分母是“该家族总预测次数”，逻辑如下：
class_stats <- class_stats %>%
  mutate(Total = Gain + Loss,
         Gain_Rate = Gain / Total,
         Loss_Rate = Loss / Total) %>%
  arrange(desc(Total)) %>% inner_join(classback,by='class')

tf_data <- class_stats %>% filter(class!='Uncharacterized')
# --- 1. 数据处理与归一化 ---
plot_data <- tf_data %>%
  # 计算全局基线 (Global Baseline Rate)
  mutate(
    global_gain_rate = sum(Gain) / sum(all),
    global_loss_rate = sum(Loss) / sum(all)
  ) %>%
  rowwise() %>%
  mutate(
    # 计算每个家族的原始发生率
    raw_rate_gain = Gain / all,
    raw_rate_loss = Loss / all,
    
    # 计算富集倍数 (O/E Ratio): 家族发生率 / 全局平均发生率
    # >1 代表该家族倾向于发生此类事件；<1 代表倾向于不发生
    Enrichment_Gain = raw_rate_gain / global_gain_rate,
    Enrichment_Loss = raw_rate_loss / global_loss_rate
  ) %>%
  ungroup() %>%
  # 转换为长格式以便绘图
  select(class, Enrichment_Gain, Enrichment_Loss) %>%
  pivot_longer(cols = c(Enrichment_Gain, Enrichment_Loss), 
               names_to = "Type", values_to = "Enrichment") %>%
  mutate(
    Type = ifelse(Type == "Enrichment_Gain", "Gain (Creating)", "Loss (Disrupting)"),
    # 按照 Gain 的富集程度对家族进行排序
    class = factor(class, levels = tf_data$class[order(tf_data$Gain / tf_data$all)])
  )

# 重新计算，把gain和loss合并
# --- 1. 数据处理与归一化 ---
plot_data <- tf_data %>%
  # 计算全局基线 (Global Baseline Rate)
  mutate(
    global_rate = sum(Total) / sum(all),
  ) %>%
  rowwise() %>%
  mutate(
    # 计算每个家族的原始发生率
    raw_rate = Total / all,

    # 计算富集倍数 (O/E Ratio): 家族发生率 / 全局平均发生率
    # >1 代表该家族倾向于发生此类事件；<1 代表倾向于不发生
    Enrichment_change = raw_rate / global_rate,
  ) %>%
  ungroup() %>%
  # 转换为长格式以便绘图
  select(class, Enrichment_change) %>%
  pivot_longer(cols = c(Enrichment_change), 
               names_to = "Type", values_to = "Enrichment") %>%
  mutate(
    Type = 'Disrupted',
    # 按照 Gain 的富集程度对家族进行排序
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
    x = "Log2 (Enrichment Ratio)", # 记得改坐标轴标签
    y = NULL
  ) +
  theme_bw()+mythem+theme(legend.position = "none")

# 1. 根据 log2FC 计算出 class 的排序顺序 (从低到高或从高到低)
# 这里我们按照 log2FC 升序排列，这样在图中就是从下往上画
class_order <- plot_data$class[order(plot_data$log2FC)]

# 2. 将两个数据框的 class 都指定为具有相同 Level 的因子
plot_data$class <- factor(plot_data$class, levels = class_order)

gc_data <- merged_data  %>% filter(class!='Uncharacterized') %>% 
  # --- 步骤 A: 计算 GC 含量 ---
  # 使用 stringr::str_count 计算 G 和 C 的数量，除以序列总长度
  mutate(gc_content = str_count(V14, "[GC]") / nchar(V14))

# 注意：这里假设 gc_data 是你的GC含量数据框
# 确保 gc_data$class 里的名称和 plot_data$class 完全一致
gc_data$class <- factor(gc_data$class, levels = class_order) 


# =======================================================
# 步骤 2: 绘制左侧棒棒糖图 (你的原代码微调)
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
  mythem + # 如果你有自定义主题请取消注释
  theme(legend.position = "none")


# =======================================================
# 步骤 3: 绘制右侧 Boxplot
# =======================================================
p2 <- ggplot(gc_data, aes(x = gc_content, y = class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "#E8F1F2") + # 箱体颜色可自选
  labs(
    x = "GC Content", 
    y = NULL
  ) +
  theme_bw() + 
  mythem + # 保持主题一致
  theme(
    # 关键：隐藏Y轴的文字、刻度和标题，因为左图已经有了
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank() # 可选：为了美观可以去掉横向网格
  )



# =======================================================
# 步骤 4: 拼接图形
# =======================================================
# 使用 patchwork 进行拼接
# widths 参数控制两个图的宽度比例，例如 c(1.5, 1) 表示左图宽一点
# 1. 修改左图 p1：去掉右边距 (r = 0)
# t=top, r=right, b=bottom, l=left
p1 <- p1 + theme(plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5, unit = "pt"))

# 2. 修改右图 p2：去掉左边距 (l = 0)
p2 <- p2 + theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0, unit = "pt"))

# 3. 拼接
combined_plot <- p1 + p2 + plot_layout(widths = c(1.5, 1.3))

# 1. 合并两个数据框，确保 class 一一对应
# inner_join 会自动根据 class 列对齐数据
correlation_data <- plot_data %>%
  select(class, log2FC) %>%
  inner_join(avg_gc_data, by = "class")
# 计算 Spearman 相关性
res <- cor.test(correlation_data$log2FC, 
                correlation_data$Mean_GC, 
                method = "spearman")

# 打印完整结果
print(res)

# 提取关键指标以便后续使用
rho_value <- round(res$estimate, 3) # 相关系数 rho
p_value   <- formatC(res$p.value, format = "e", digits = 2) # P值 (科学计数法)

# 输出一句话结果
paste0("Spearman Rho = ", rho_value, ", P-value = ", p_value)

p_corr <- ggplot(correlation_data, aes(x = log2FC, y = Mean_GC)) +
  geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
  geom_smooth(method = "lm", color = "red", se = FALSE, linetype = "dashed") + # 添加拟合线
  
  # 自动添加相关系数和P值标签
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  
  labs(
    x = "Log2 (Enrichment Ratio)",
    y = "Average GC Content",
  ) +
  theme_bw()+mythem
avg_gc_data <- gc_data %>%
  group_by(class) %>%
  summarise(Mean_GC = mean(gc_content, na.rm = TRUE)) # 计算均值

# 关键：确保右图数据的 class 也是同样的因子顺序
avg_gc_data$class <- factor(avg_gc_data$class, levels = class_order)
p2 <- ggplot(avg_gc_data, aes(x = Mean_GC, y = class)) +
  geom_col(width = 0.6, fill = "#7F8C8D", alpha = 0.8) + # 使用中性灰色，宽度调细一点更精致
  
  # 如果想在柱子旁边显示具体数值，可以加上下面这行：
  # geom_text(aes(label = round(Mean_GC, 2)), hjust = -0.2, size = 3) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) + # 让柱子紧贴Y轴，右边留点空隙
  labs(x = "Average GC Content", y = NULL) +
  theme_bw() +
  theme(
    # 关键：隐藏 Y 轴所有元素，实现无缝拼接
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(), # 去掉横向网格让画面更干净
    plot.margin = margin(l = 0) # 去掉左边距，让两个图靠得更近
  )


# =======================================================
# 步骤 4: 拼接
# =======================================================
# 调整 widths 比例，通常 Enrichment 图为主，GC 图为辅，可以设为 2:1
combined_plot <- p1 + p2 + plot_layout(widths = c(1.8, 1))
# --- 2. 绘制棒棒糖图 (Lollipop Chart) ---
ggplot(plot_data, aes(x = Enrichment, y = class, color = Type)) +
  # 添加基准线 x=1 (代表平均水平)
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  
  # 绘制连接线
  geom_line(aes(group = class), color = "grey80", size = 0.5) +
  
  # 绘制点
  geom_point(size = 3.5) +
  
  # 配色 (Gain用红色系，Loss用蓝色系，符合常规直觉)
  scale_color_manual(values = c("Gain (Creating)" = "#D6604D", "Loss (Disrupting)" = "#4393C3")) +
  
  labs(
    title = "Relative Enrichment of MNV-induced TFBS Changes",
    subtitle = "Normalized by total motif count (O/E Ratio)",
    x = "Fold Enrichment (Observed Rate / Global Average)",
    y = NULL,
    color = "Event Type"
  ) +
  
  theme_bw() +mythem

# 准备绘图数据（长格式）
plot_data <- class_stats %>%
  select(class, Gain, Loss) %>%
  pivot_longer(cols = c(Gain, Loss), names_to = "Effect", values_to = "Count")



# 绘制堆叠柱状图显示不同 TF Class 的受影响程度
ggplot(plot_data, aes(x = reorder(class, Count), y = Count, fill = Effect)) +
  geom_bar(stat = "identity") +
  coord_flip() + # 横向显示，方便阅读较长的家族名称
  scale_fill_manual(values = c("Gain" = "#D6604D", "Loss" = "#438fa9")) +
  labs(title = "Impact of MNVs on different TF Classes",
       x = "TF Class (JASPAR)",
       y = "Number of affected TFBS") +
  theme_bw()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          