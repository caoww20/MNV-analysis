# 功能：评估心肌病 GWAS 信号在高 FST MNV 区域内的富集，使用超几何/二项/泊松检验计算显著性。
# 用法：直接运行，依赖输入路径已写死；如需复用请替换 FST 列表、MNV 坐标及 GWAS 位点文件。
a <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/mnvid_usefor_mapping_high_fst',col_names = F)
b <- read_tsv('/home/caow/03mnv/analyse3/02pop/02group/03fst/significant/mnv_sig.txt')
colnames(a)[1:2] <- colnames(b)[1:2]
fstmnv <- b %>% inner_join(a,by=c('CHROM','POS')) %>% distinct(X3)
colnames(fstmnv) <- 'mnvid'
mnvpos1000g <- read_tsv('/home/caow/03mnv/analyse3/00addition/data/1000gmnv.pos.bed',col_names = F)
mnvpos1000g[[1]] <- sub("^chr", "", mnvpos1000g[[1]])

colnames(mnvpos1000g)[c(1,3,4)] <- c('CHROM','POS','mnvid')
fstmnvpos <- fstmnv %>% inner_join(mnvpos1000g,by='mnvid') %>% mutate(CHROM=as.double(CHROM))

cardi_gwas <- read_tsv('/home/caow/03mnv/analyse3/02pop/02group/05gwas/dilated_cardiomyopathy.gwas',col_names = F) %>% select(c(12,13))
colnames(cardi_gwas) <- c('CHROM','POS')
cardi_gwas %>% inner_join(fstmnvpos,by=c('CHROM','POS'))
                          

allsnv <- read_tsv('/home/caow/03mnv/analyse3/02data/03density/02datasets/snv_1000G.txt',col_names = F)
colnames(allsnv)[2:3] <- c('CHROM','POS')
allsnv %>% inner_join(fstmnvpos,by=c('CHROM','POS'))


n_total_snv <- 106     # 所有的心肌病GWAS相关SNV数
n_observed  <- 12        # 心肌病GWAS相关SNV落在高FST的MNV区域内的SNV数

total_snvs_background <- 61599150  # 背景总SNV 
target_snvs_background <- 76892    # 背景中落在高FST的MNV区域的SNV

background_prob <- target_snvs_background / total_snvs_background

# 2. 计算期望值 lambda
lambda <- n_total_snv * background_prob
print(paste("Expected number (Lambda):", lambda))

# 3. 计算泊松分布概率 P(X >= 12)
# ppois 计算的是 P(X <= k)，所以我们要用 1 - ppois(k-1) 或者 lower.tail = FALSE
p_value <- ppois(n_observed - 1, lambda = lambda, lower.tail = FALSE)

print(paste("P-value:", p_value))


# P(X >= 12 | n=100, p=background_prob)
binom.test(x = 12, n = 106, p = background_prob, alternative = "greater")

phyper(11, 76892, 61599150 - 76892, 106, lower.tail = FALSE)
