# Purpose: explore MNV effects on TFBS and target-gene constraint, compute nearest TSS,
# join gnomAD LOEUF/pLI, and plot distributions and proportions.
# Usage: source after preparing tfbs_res, jaspar_meta, and gnomad constraint tables;
# assumes hg38 GFF3 at ~/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3.
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer) # For reading bigWig (Layer 3)

# 1. Load gnomAD gene constraint data (download required)
# File assumed to be gnomad_constraint.txt
# Key columns: gene, oe_lof_upper (LOEUF), pLI
gnomad <- read_tsv("/home/caow/03mnv/analyse3/00addition/data/gnomad.v2.1.1.lof_metrics.by_gene.txt") %>%
  select(gene, LOEUF = oe_lof_upper, pLI)

# --- Mock gnomAD data (for demo only) ---
# set.seed(123)
# all_genes <- unique(c(tf_data$class, "TP53", "MYC", "GAPDH")) # Example gene names
# gnomad <- data.frame(
#   gene = paste0("Gene", 1:1000),
#   LOEUF = runif(1000, 0, 2), # 0-0.35 is constrained
#   pLI = runif(1000, 0, 1)    # >0.9 is constrained
# )
# # ----------------------------------------

# 2. Load your MNV-TFBS results (jaspar_pair.txt)
# Assumed columns: MNVID, TF_Matrix_ID, Effect(Gain/Loss), MNV_Chr, MNV_Pos
# mnv_tfbs <- read.table("jaspar_pair.txt", header=FALSE)
# colnames(mnv_tfbs)[c(1, 10)] <- c("matrix_id", "effect")
# MNV coordinates are needed; assumed in cols 4-5 or join back
mnv_tfbs <- tfbs_res
mnv_tfbs[[3]] <- sub("^chr", "", mnv_tfbs[[3]])

# 3. Load JASPAR metadata (map Matrix ID to TF Gene Name)
# jaspar_meta <- read.table("jaspar_metadata.tsv", header=TRUE, sep="\t", quote="") %>%
  # select(matrix_id, name, symbol = uniprot_ids) # Convert ID to Gene Symbol if needed

jaspar_meta <- jaspar_meta %>% select(matrix_id, name, symbol = uniprot_ids)
# --- Step 1: define Target Gene (nearest TSS) ---

# A. Prepare GRanges for TFBS-MNV
# Assumes mnv_tfbs has chr, start, end
gr_mnv <- GRanges(seqnames = mnv_tfbs$V3, 
                  ranges = IRanges(start = mnv_tfbs$V12, end = mnv_tfbs$V13))

# B. Prepare GRanges for gene TSS (from GFF3)
gff <- import.gff3("~/01miRNASNP/gene/Homo_sapiens.GRCh38.110.gff3")
genes <- gff[gff$type == "gene"]
promoters <- promoters(genes, upstream=0, downstream=1) # Get TSS
gr_tss <- promoters

# --- Mock TSS data ---
# gr_tss <- GRanges(seqnames = sample(paste0("chr", 1:22), 1000, replace=T),
#                   ranges = IRanges(start = sample(1:1e8, 1000), width=1),
#                   gene_name = paste0("Gene", 1:1000))

# C. Find nearest TSS
nearest_idx <- distanceToNearest(gr_mnv, gr_tss)
mnv_tfbs$target_gene <- mcols(gr_tss)$Name[subjectHits(nearest_idx)]
mnv_tfbs$dist_to_tss <- mcols(nearest_idx)$distance

# --- Step 2: join constraint scores ---

analysis_layer1 <- mnv_tfbs %>%
  left_join(gnomad, by = c("target_gene" = "gene")) %>%
  filter(!is.na(LOEUF)) %>% # Drop genes without scores
  mutate(
    # Define constraint category
    constraint_bin = case_when(
      LOEUF <= 1.48 ~ "Constrained",
      LOEUF > 1.48 ~ "Tolerant"
    )
  )

# --- Step 3: statistics and plots (key) ---
# Compare: fraction of constrained genes among MNV-disrupted genes vs background
# Or: observed / expected

# Quick visualization: LOEUF distribution for disrupted TFBS targets vs background
ggplot(analysis_layer1, aes(x = LOEUF)) +
  geom_density(aes(fill = "MNV Targets"), alpha = 0.4) +
  geom_density(data = gnomad, aes(x = LOEUF, fill = "Genome Background"), alpha = 0.4) +
  labs(title = "Layer 1: MNV Targets are skewed towards tolerant genes",
       x = "LOEUF Score (Lower = More Constrained)") +
  theme_bw()

# Statistical test
wilcox.test(analysis_layer1$LOEUF, gnomad$LOEUF, alternative = "greater")
# If p < 0.05 and MNV LOEUF mean is higher (less constrained),
# MNVs tend to occur near less essential genes.

library(dplyr)
library(tidyr)

# --- A. Load full MNV table (your raw MNV file) ---
# Must include: mnvid, target_gene
# Assume object exists or read from file
# all_mnvs <- read.table("mnv_in_utr5.txt", ...) 
all_mnvs <- read_tsv('~/mnv2tfbs/hg38mnv_promoter',col_names = F,col_types = cols(
  X2 = col_character(),        # positions: must be character (contains commas)
)) %>%
  transmute(
    chr = X1,
    start = as.integer(str_extract(X2, "^[0-9]+")),   # First number
    end   = as.integer(str_extract(X2, "[0-9]+$")),    # Last number
    mnvid = X3
  )
all_mnv_r <- GRanges(seqnames = all_mnvs$chr, 
                  ranges = IRanges(start = all_mnvs$start, end = all_mnvs$end))
nearest_idx <- distanceToNearest(all_mnv_r, gr_tss)
all_mnvs$target_gene <- mcols(gr_tss)$Name[subjectHits(nearest_idx)]
all_mnvs$dist_to_tss <- mcols(nearest_idx)$distance

# --- B. Load TFBS effects table (your TFBS analysis) ---
# Must include: mnvid, effect (Gain/Loss)
# tfbs_effects <- read.table("jaspar_analysis_result.txt", ...)
tfbs_effects <- mnv_tfbs %>% select(V5,effect)
colnames(tfbs_effects)[1] <- 'mnvid'  

mnv_final_data <- all_mnvs %>%
  # 1. Join TFBS results
  left_join(tfbs_effects, by = "mnvid") %>%
  
  # 2. Group by MNV to collapse records
  group_by(mnvid) %>%
  summarise(
    # Keep target_gene (assume one per MNV; take first)
    target_gene = target_gene[1],
    
    # 3. Priority logic
    # Check all effect records per MNV
    category = case_when(
      # Priority 1: if any "Loss", classify as "Disrupting (Loss)"
      # (even if it also has Gain; disruption is more severe)
      any(effect == "Loss", na.rm = TRUE) ~ "Disrupting (Loss)",
      
      # Priority 2: if no Loss but has Gain, classify as "Creating (Gain)"
      any(effect == "Gain", na.rm = TRUE) ~ "Creating (Gain)",
      
      # Priority 3: neither Loss nor Gain (NA), classify as "Neutral"
      TRUE ~ "Neutral (Background)"
    ),
    .groups = "drop" # Ungroup
  )

# --- D. Join gnomAD constraint data ---
# Ensure gnomad is loaded and thresholds computed
loeuf_threshold <- quantile(gnomad$LOEUF, probs = 0.20, na.rm = TRUE)
constrained_genes <- gnomad$gene[gnomad$LOEUF <= loeuf_threshold]

plot_dataset <- mnv_final_data %>%
  filter(target_gene %in% gnomad$gene) %>% # Drop missing LOEUF data
  mutate(
    in_constrained = target_gene %in% constrained_genes
  )

# --- E. Define calculation helper (with guards) ---
calc_fraction <- function(df, cat_name) {
  n_total <- nrow(df)
  n_cons <- sum(df$in_constrained, na.rm=TRUE)
  
  if(n_total == 0) return(NULL) # Skip empty group
  
  ci <- binom.test(n_cons, n_total)$conf.int
  data.frame(
    Category = cat_name,
    Fraction = n_cons / n_total,
    Lower = ci[1], 
    Upper = ci[2],
    n_size = n_total
  )
}

# --- F. Build plot data ---
final_plot_data <- rbind(
  calc_fraction(plot_dataset %>% filter(category == "Disrupting (Loss)"), "Disrupting (Loss)"),
  calc_fraction(plot_dataset %>% filter(category == "Creating (Gain)"), "Creating (Gain)"),
  calc_fraction(plot_dataset %>% filter(category == "Neutral (Background)"), "Unchaged")
)

# --- G. Plot (Wang et al. style) ---
library(ggplot2)

ggplot(final_plot_data, aes(x = Category, y = Fraction, color = Category)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, size = 1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 0.3)) +
  
  # Colors: Loss (red/orange), Gain (green/blue), Neutral (gray)
  scale_color_manual(values = c(
    "Disrupting (Loss)" = "#D55E00",       # Highlight
    "Creating (Gain)" = "#009E73",         # Contrast
    "Unchanged" = "grey50"      # Background
  )) +
  
  labs(
    y = "Fraction of MNVs on Constrained Genes",
    x = NULL
  ) +
  theme_bw() +
  mythem+theme(
    legend.position = "none",
  )


# --- Step 1: join TF gene constraint ---

# Note: jaspar_meta uniprot_ids may need conversion to Gene Symbol for gnomad
# Assume jaspar_meta already has symbol column

tf_constraint <- jaspar_meta %>%
  # Split motifs that map to multiple genes
  separate_rows(name, sep = "::") %>% 
  mutate(name = toupper(name)) %>% 
  left_join(gnomad %>% distinct(), by = c("name" = "gene")) %>%
  group_by(matrix_id) %>%
  # If a motif maps to multiple TFs, use the most constrained (min LOEUF)
  summarise(tf_LOEUF = min(LOEUF, na.rm = TRUE))
