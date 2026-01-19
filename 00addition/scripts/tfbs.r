# 功能：可视化 TFBS 相关的期望/长度分布、GC 含量、重复序列比例等；依赖前置数据框 df、expected_var、merged_data 等已在环境中加载。
# 用法：在准备好数据（见脚本内变量命名）后 source 本脚本，生成各类 ggplot 可视化。
expvar <- df %>% left_join(expected_var,by='type') %>% select(c(1,8))
expvar[,2] <- expvar[,2]/mnvRegion[,2]
expvar <- expvar[-15,]
expvar$id<-1:nrow(expvar)
res <- reshape2::melt(expvar,id=c('type','id')) 

ggplot(res, aes(x = id, y = value, color = variable)) +
  geom_point(size=2) +
  geom_line(lwd=1) +
  labs(title='',x=NULL,y='Expected variants/Length')+
  scale_x_continuous(breaks=c(1:nrow(df_ratio)),labels=df_ratio$type)+
  theme_bw()+mythem+theme(axis.text.x = element_text(angle = 45,hjust = 1))

# TFBS 计算 GC 含量 & 绘图
# =======================================================

merged_data  %>% filter(class!='Uncharacterized') %>% 
  # --- 步骤 A: 计算 GC 含量 ---
  # 使用 stringr::str_count 计算 G 和 C 的数量，除以序列总长度
  mutate(gc_content = str_count(V14, "[GC]") / nchar(V14)) %>%
  
  # --- 步骤 B: 绘图 ---
  # reorder(class, gc_content): 让箱线图按 GC 含量的中位数自动排序，看起来更整洁
  ggplot(aes(x = reorder(class, gc_content), y = gc_content, fill = class)) +
  
  # 绘制箱线图
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + # alpha设置透明度
  
  # # 添加抖动点 (Jitter): 如果样本量少，建议加上这个，能看到每个点及其分布
  # geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  
  # 翻转坐标轴 (因为分类名称太长，放在 Y 轴更好看)
  coord_flip() +
  
  # 设置标签
  labs(
    title = "GC Content Distribution by TF Class",
    x = "TF Class",
    y = "GC Content"
  ) +
  
  # 美化主题
  theme_bw() +mythem+
  theme(
    legend.position = "none", # 隐藏图例（因为Y轴已经是分类了）
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

library(stringr) # 需要先安装 install.packages("stringr")


# seq: 输入的序列 (V14)
# start: 序列起始位置 (V12)
# pos_str: 突变绝对位置 (V4)
# allele_str: 要替换进去的碱基 (如果是Loss用V7，如果是Gain用V6)
replace_allele <- function(seq, start, pos_str, allele_str) {
  
  # 将序列拆分为字符向量
  seq_chars <- str_split(seq, "")[[1]]
  
  # 解析位置和碱基
  positions <- as.numeric(unlist(str_split(pos_str, ",")))
  alleles   <- unlist(str_split(allele_str, ","))
  
  # 遍历替换
  for (i in seq_along(positions)) {
    pos <- positions[i]
    target_allele <- alleles[i]
    
    # 计算相对位置
    rel_idx <- pos - start + 1
    
    # 安全检查并替换
    if (rel_idx > 0 && rel_idx <= length(seq_chars)) {
      seq_chars[rel_idx] <- target_allele
    }
  }
  
  # 重新组合
  return(paste(seq_chars, collapse = ""))
}

# B. 重复区域判断函数 (将您的逻辑封装在这里)
check_is_repeat <- function(seq_vector) {
  # 逻辑 A: 均聚物 (5个及以上相同碱基)
  is_homopolymer <- str_detect(seq_vector, "([A-Z])\\1{4,}")
  
  # 逻辑 B: 双碱基重复 (4个及以上单元，如 ACACACAC)
  is_di_repeat <- str_detect(seq_vector, "([A-Z]{2})\\1{3,}")
  
  # 返回综合结果 (只要满足其中一个即为 TRUE)
  return(is_homopolymer | is_di_repeat)
}

plot_data <- merged_data%>% filter(effect=='Loss') %>%
  rowwise() %>%
  mutate(
    # --- 第一步：生成序列 (Ref 和 Alt) ---
    Ref_Seq_Final = ifelse(effect == "Loss", V14, replace_allele(V14, V12, V4, V6)),
    Alt_Seq_Final = ifelse(effect == "Gain", V14, replace_allele(V14, V12, V4, V7))
  ) %>%
  ungroup() %>%
  
  # --- 第二步：同时应用重复判断逻辑 ---
  # 这里不需要 rowwise，因为 check_is_repeat 函数支持向量化操作，速度更快
  mutate(
    # 判断 Ref 序列是否由重复构成
    Ref_is_Repeat = check_is_repeat(Ref_Seq_Final),
    
    # 判断 Alt 序列是否由重复构成
    Alt_is_Repeat = check_is_repeat(Alt_Seq_Final),
    
    # --- 第三步：综合标记 ---
    # 通常如果 Ref 或 Alt 任意一个是重复区域，该位点都可能被视为"低复杂度区域" (Low Complexity)
    Is_Low_Complexity = Ref_is_Repeat | Alt_is_Repeat
  ) %>% group_by(class) %>%
  summarise(
    total_count = n(),                          # 该家族总位点数
    repeat_count = sum(Is_Low_Complexity),      # 重复位点数 (TRUE的数量)
    repeat_rate = mean(Is_Low_Complexity)       # 重复率 (TRUE的比例)
  ) %>%
  ungroup()

ggplot(plot_data, aes(x = reorder(class, repeat_rate), y = repeat_rate)) +
  
  # A. 绘制柱子
  geom_col(aes(fill = repeat_rate), width = 0.7) +
  
  # B. 添加文字标签 (百分比 + 样本量)
  # label 格式示例: "75.0% (n=20)"
  geom_text(aes(label = paste0(sprintf("%.1f", repeat_rate * 100), "%", 
                               "\n(n=", total_count, ")")), 
            hjust = -0.1,  # 文字放在柱子右侧
            size = 3.5, 
            color = "black") +
  
  # C. 坐标轴翻转 (名字太长时必备)
  coord_flip() +
  
  # D. 设置Y轴范围 (给文字标签留出空间，防止溢出)
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0, max(plot_data$repeat_rate) * 1.2)) +
  
  # E. 颜色渐变 (可选，颜色越深代表重复率越高)
  scale_fill_gradient(low = "#85D4E3", high = "#F4B5BD") +
  
  # F. 标签与主题
  labs(
    title = "Proportion of Repetitive/Low-Complexity Sites by TF Family",
    x = "",  # Y轴已经很清楚是分类名了，X轴标签可以留空
    y = "Percentage of Repetitive Sequences"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",              # 隐藏图例
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5) # 标题居中
  )
