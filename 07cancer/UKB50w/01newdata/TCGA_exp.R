mythem=theme(panel.grid=element_blank(),
             plot.margin=unit(rep(2,4),'lines'),
             legend.position="none",
             panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
             # linewidth replaces deprecated element_rect
             # title = element_text(vjust=1,size=16,face = "bold", color='black'),
             axis.text.x=element_text(vjust=1,size=8,face = "bold", colour='black'),
             axis.title.x=element_text(vjust=0, size=8,face = "bold", color='black'),
             axis.text.y=element_text(vjust=0.5,size=8,face = "bold", color = 'black'),
             axis.title.y=element_text(vjust=2, size=8,face = "bold", color='black')
)

# setwd('/data/jinww/04reference/publicDB/PancanMNVQTL/phenotype/eQTL')
# LUAD=read.table('LUAD_eQTL_quantity.filtered',header = T,sep = '\t')
setwd("/data/jinww/04reference/publicDB/TCGA/expression_20220422")
LUAD2=read.table('LUAD_gene_exp',header = T,sep = '\t')
df1 <- LUAD2[, grep("\\.01", names(LUAD2))]
df2 <- LUAD2[, grep("\\.11", names(LUAD2))]


cancer_exp=t(LUAD2[LUAD2$gene=='ENSG00000164741.15|DLC1',grep("\\.01", names(LUAD2))])
normal_exp=t(LUAD2[LUAD2$gene=='ENSG00000164741.15|DLC1',grep("\\.11", names(LUAD2))])

cancer_exp=t(LUAD2[LUAD2$gene=='ENSG00000141338.14|ABCA8',grep("\\.01", names(LUAD2))])
normal_exp=t(LUAD2[LUAD2$gene=='ENSG00000141338.14|ABCA8',grep("\\.11", names(LUAD2))])

data <- data.frame(
  group = c(rep(c("cancer"),nrow(cancer_exp)),rep(c("normal"),nrow(normal_exp))),
  value = c(cancer_exp[,1], normal_exp[,1])
)


p=ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  labs(title = "", x = "", y = "The expression of DLC1 (FPKM)")+
  theme_bw()+mythem
ggsave(
  filename = 'ABCA8.pdf',
  plot = p,width = 4,height = 5,
  units = 'in'
)


