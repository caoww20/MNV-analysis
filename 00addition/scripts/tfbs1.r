# 功能：探索 MNV 对 TFBS 及靶基因约束力的影响，计算最近 TSS、关联 gnomAD LOEUF/pLI，并绘制分布与比例。
# 用法：在准备好 tfbs_res、jaspar_meta、gnomad 约束表等对象后 source 本脚本；脚本假设 hg38 GFF3 位于 ~/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3。
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer) # 如果需要读取 bigWig (Layer 3)

# 1. 读取 gnomAD 基因约束力数据 (你需要下载这个文件)
# 假设文件名为 gnomad_constraint.txt
# 重点列: gene, oe_lof_upper (即 LOEUF), pLI
gnomad <- read_tsv("/home/caow/03mnv/analyse3/00addition/data/gnomad.v2.1.1.lof_metrics.by_gene.txt") %>%
  select(gene, LOEUF = oe_lof_upper, pLI)

# --- 模拟 gnomAD 数据 (仅供代码跑通演示) ---
# set.seed(123)
# all_genes <- unique(c(tf_data$class, "TP53", "MYC", "GAPDH")) # 假设这是基因名
# gnomad <- data.frame(
#   gene = paste0("Gene", 1:1000),
#   LOEUF = runif(1000, 0, 2), # 0-0.35 为 constrained
#   pLI = runif(1000, 0, 1)    # >0.9 为 constrained
# )
# # ----------------------------------------

# 2. 读取你的 MNV-TFBS 结果 (jaspar_pair.txt)
# 假设结构: MNVID, TF_Matrix_ID, Effect(Gain/Loss), MNV_Chr, MNV_Pos
# mnv_tfbs <- read.table("jaspar_pair.txt", header=FALSE)
# colnames(mnv_tfbs)[c(1, 10)] <- c("matrix_id", "effect")
# 还需要 MNV 的坐标，假设在第 4-5 列或者你需要关联回去
mnv_tfbs <- tfbs_res
mnv_tfbs[[3]] <- sub("^chr", "", mnv_tfbs[[3]])

# 3. 读取 JASPAR 元数据 (关联 Matrix ID 到 TF Gene Name)
# jaspar_meta <- read.table("jaspar_metadata.tsv", header=TRUE, sep="\t", quote="") %>%
  # select(matrix_id, name, symbol = uniprot_ids) # 注意: 这里需要把 ID 转为 Gene Symbol

jaspar_meta <- jaspar_meta %>% select(matrix_id, name, symbol = uniprot_ids)
# --- 步骤 1: 定义 Target Gene (Nearest TSS) ---

# A. 准备 TFBS-MNV 的 GRanges 对象
# 假设 mnv_tfbs 有 chr, start, end
gr_mnv <- GRanges(seqnames = mnv_tfbs$V3, 
                  ranges = IRanges(start = mnv_tfbs$V12, end = mnv_tfbs$V13))

# B. 准备基因 TSS 的 GRanges 对象 (从 GFF3 提取)
gff <- import.gff3("~/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3")
genes <- gff[gff$type == "gene"]
promoters <- promoters(genes, upstream=0, downstream=1) # 获取 TSS
gr_tss <- promoters

# --- 模拟 TSS 数据 ---
# gr_tss <- GRanges(seqnames = sample(paste0("chr", 1:22), 1000, replace=T),
#                   ranges = IRanges(start = sample(1:1e8, 1000), width=1),
#                   gene_name = paste0("Gene", 1:1000))

# C. 寻找最近的 TSS
nearest_idx <- distanceToNearest(gr_mnv, gr_tss)
mnv_tfbs$target_gene <- mcols(gr_tss)$Name[subjectHits(nearest_idx)]
mnv_tfbs$dist_to_tss <- mcols(nearest_idx)$distance

# --- 步骤 2: 关联约束力分数 ---

analysis_layer1 <- mnv_tfbs %>%
  left_join(gnomad, by = c("target_gene" = "gene")) %>%
  filter(!is.na(LOEUF)) %>% # 过滤掉没有分数的基因
  mutate(
    # 定义约束类别
    constraint_bin = case_when(
      LOEUF <= 1.48 ~ "Constrained",
      LOEUF > 1.48 ~ "Tolerant"
    )
  )

# --- 步骤 3: 统计与绘图 (关键) ---
# 我们比较: 发生 MNV disruption 的基因中，Constrained 的比例 vs 基因组背景比例
# 或者: 比较 Observed / Expected

# 简单可视化: 比较 Disrupted TFBS Target Genes 的 LOEUF 分布 vs 所有基因背景
ggplot(analysis_layer1, aes(x = LOEUF)) +
  geom_density(aes(fill = "MNV Targets"), alpha = 0.4) +
  geom_density(data = gnomad, aes(x = LOEUF, fill = "Genome Background"), alpha = 0.4) +
  labs(title = "Layer 1: MNV Targets are skewed towards tolerant genes",
       x = "LOEUF Score (Lower = More Constrained)") +
  theme_bw()

# 统计检验
wilcox.test(analysis_layer1$LOEUF, gnomad$LOEUF, alternative = "greater")
# 如果 p < 0.05 且 MNV 组 LOEUF 均值更高（更不保守），说明 MNV 倾向于发生在不重要的基因旁（即重要基因旁的被清除了）。

library(dplyr)
library(tidyr)

# --- A. 模拟读取全量 MNV 表 (你的原始 MNV 文件) ---
# 必须包含: mnvid, target_gene
# 假设你已经有了这个对象，或者从文件读取
# all_mnvs <- read.table("mnv_in_utr5.txt", ...) 
all_mnvs <- read_tsv('~/mnv2tfbs/hg38mnv_promoter',col_names = F,col_types = cols(
  X2 = col_character(),        # positions：必须是字符（含逗号）
)) %>%
  transmute(
    chr = X1,
    start = as.integer(str_extract(X2, "^[0-9]+")),   # 第一个数字
    end   = as.integer(str_extract(X2, "[0-9]+$")),    # 最后一个数字
    mnvid = X3
  )
all_mnv_r <- GRanges(seqnames = all_mnvs$chr, 
                  ranges = IRanges(start = all_mnvs$start, end = all_mnvs$end))
nearest_idx <- distanceToNearest(all_mnv_r, gr_tss)
all_mnvs$target_gene <- mcols(gr_tss)$Name[subjectHits(nearest_idx)]
all_mnvs$dist_to_tss <- mcols(nearest_idx)$distance

# --- B. 模拟读取 TFBS 效应表 (你的 TFBS 分析结果) ---
# 必须包含: mnvid, effect (Gain/Loss)
# tfbs_effects <- read.table("jaspar_analysis_result.txt", ...)
tfbs_effects <- mnv_tfbs %>% select(V5,effect)
colnames(tfbs_effects)[1] <- 'mnvid'  

mnv_final_data <- all_mnvs %>%
  # 1. 连接 TFBS 结果
  left_join(tfbs_effects, by = "mnvid") %>%
  
  # 2. 按 MNV 分组进行“降维”聚合
  group_by(mnvid) %>%
  summarise(
    # 保留 target_gene (假设每个 MNV 只有一个靶基因，取第一个即可)
    target_gene = target_gene[1],
    
    # 3. 【核心修改】优先级判定逻辑
    # 检查该 MNV 下的所有 effect 记录
    category = case_when(
      # 优先级 1: 只要包含 "Loss"，就定性为 "Disrupting (Loss)"
      # (即使它同时也包含 Gain，我们认为破坏原有功能的后果更严重)
      any(effect == "Loss", na.rm = TRUE) ~ "Disrupting (Loss)",
      
      # 优先级 2: 如果没有 Loss 但包含 "Gain"，定性为 "Creating (Gain)"
      any(effect == "Gain", na.rm = TRUE) ~ "Creating (Gain)",
      
      # 优先级 3: 既没 Loss 也没 Gain (即 NA)，定性为 "Neutral"
      TRUE ~ "Neutral (Background)"
    ),
    .groups = "drop" # 解除分组
  )

# --- D. 关联 gnomAD 约束力数据 ---
# 确保 gnomad 数据已加载且计算了阈值
loeuf_threshold <- quantile(gnomad$LOEUF, probs = 0.20, na.rm = TRUE)
constrained_genes <- gnomad$gene[gnomad$LOEUF <= loeuf_threshold]

plot_dataset <- mnv_final_data %>%
  filter(target_gene %in% gnomad$gene) %>% # 过滤掉没有 LOEUF 数据的
  mutate(
    in_constrained = target_gene %in% constrained_genes
  )

# --- E. 定义计算函数 (带容错) ---
calc_fraction <- function(df, cat_name) {
  n_total <- nrow(df)
  n_cons <- sum(df$in_constrained, na.rm=TRUE)
  
  if(n_total == 0) return(NULL) # 跳过空组
  
  ci <- binom.test(n_cons, n_total)$conf.int
  data.frame(
    Category = cat_name,
    Fraction = n_cons / n_total,
    Lower = ci[1], 
    Upper = ci[2],
    n_size = n_total
  )
}

# --- F. 生成绘图数据 ---
final_plot_data <- rbind(
  calc_fraction(plot_dataset %>% filter(category == "Disrupting (Loss)"), "Disrupting (Loss)"),
  calc_fraction(plot_dataset %>% filter(category == "Creating (Gain)"), "Creating (Gain)"),
  calc_fraction(plot_dataset %>% filter(category == "Neutral (Background)"), "Unchaged")
)

# --- G. 绘图 (仿 Wang et al.) ---
library(ggplot2)

ggplot(final_plot_data, aes(x = Category, y = Fraction, color = Category)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, size = 1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 0.3)) +
  
  # 颜色: Loss(红/橙), Gain(绿/蓝), Neutral(灰)
  scale_color_manual(values = c(
    "Disrupting (Loss)" = "#D55E00",       # 醒目色
    "Creating (Gain)" = "#009E73",         # 对比色
    "Unchanged" = "grey50"      # 背景色
  )) +
  
  labs(
    y = "Fraction of MNVs on Constrained Genes",
    x = NULL
  ) +
  theme_bw() +
  mythem+theme(
    legend.position = "none",
  )


# --- 步骤 1: 关联 TF 基因的约束力 ---

# 注意: jaspar_meta 的 uniprot_ids 可能需要转为 Gene Symbol 才能和 gnomad 对应
# 这里假设 jaspar_meta 已经有了 symbol 列

tf_constraint <- jaspar_meta %>%
  # 拆分可能对应多个基因的 motif
  separate_rows(name, sep = "::") %>% 
  mutate(name = toupper(name)) %>% 
  left_join(gnomad %>% distinct(), by = c("name" = "gene")) %>%
  group_by(matrix_id) %>%
  # 如果一个 motif 对应多个 TF，取最保守的那个 (min LOEUF) 作为该 motif 的约束力
  summarise(tf_LOEUF = min(LOEUF, na.rm = TRUE))
