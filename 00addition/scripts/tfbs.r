# Purpose: visualize TFBS expected/length distributions, GC content, repeat ratios, etc.; requires df, expected_var, merged_data preloaded in the environment.
# Usage: source this script after preparing data (see variable names) to generate ggplot visualizations.
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

# TFBS GC content calculation and plotting
# =======================================================

merged_data  %>% filter(class!='Uncharacterized') %>% 
  # --- Step A: compute GC content ---
  # Use stringr::str_count for G and C, divide by sequence length
  mutate(gc_content = str_count(V14, "[GC]") / nchar(V14)) %>%
  
  # --- Step B: plot ---
  # reorder(class, gc_content): sort boxplots by median GC content
  ggplot(aes(x = reorder(class, gc_content), y = gc_content, fill = class)) +
  
  # Draw boxplot
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + # alpha controls transparency
  
  # # Add jitter points: useful when sample size is small
  # geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  
  # Flip coordinates (long category labels)
  coord_flip() +
  
  # Labels
  labs(
    title = "GC Content Distribution by TF Class",
    x = "TF Class",
    y = "GC Content"
  ) +
  
  # Theme tweaks
  theme_bw() +mythem+
  theme(
    legend.position = "none", # Hide legend (categories already on axis)
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

library(stringr) # Install first via install.packages("stringr")


# seq: input sequence (V14)
# start: sequence start (V12)
# pos_str: absolute mutation positions (V4)
# allele_str: replacement alleles (V7 for Loss, V6 for Gain)
replace_allele <- function(seq, start, pos_str, allele_str) {
  
  # Split sequence into character vector
  seq_chars <- str_split(seq, "")[[1]]
  
  # Parse positions and alleles
  positions <- as.numeric(unlist(str_split(pos_str, ",")))
  alleles   <- unlist(str_split(allele_str, ","))
  
  # Iterate replacements
  for (i in seq_along(positions)) {
    pos <- positions[i]
    target_allele <- alleles[i]
    
    # Compute relative position
    rel_idx <- pos - start + 1
    
    # Bounds check and replace
    if (rel_idx > 0 && rel_idx <= length(seq_chars)) {
      seq_chars[rel_idx] <- target_allele
    }
  }
  
  # Recombine
  return(paste(seq_chars, collapse = ""))
}

# B. Repeat-region check function (encapsulates logic)
check_is_repeat <- function(seq_vector) {
  # Logic A: homopolymer (5+ identical bases)
  is_homopolymer <- str_detect(seq_vector, "([A-Z])\\1{4,}")
  
  # Logic B: dinucleotide repeat (4+ units, e.g. ACACACAC)
  is_di_repeat <- str_detect(seq_vector, "([A-Z]{2})\\1{3,}")
  
  # Return combined result (TRUE if either condition holds)
  return(is_homopolymer | is_di_repeat)
}

plot_data <- merged_data%>% filter(effect=='Loss') %>%
  rowwise() %>%
  mutate(
    # --- Step 1: build sequences (Ref and Alt) ---
    Ref_Seq_Final = ifelse(effect == "Loss", V14, replace_allele(V14, V12, V4, V6)),
    Alt_Seq_Final = ifelse(effect == "Gain", V14, replace_allele(V14, V12, V4, V7))
  ) %>%
  ungroup() %>%
  
  # --- Step 2: apply repeat logic ---
  # rowwise not needed; check_is_repeat is vectorized
  mutate(
    # Whether Ref sequence is repetitive
    Ref_is_Repeat = check_is_repeat(Ref_Seq_Final),
    
    # Whether Alt sequence is repetitive
    Alt_is_Repeat = check_is_repeat(Alt_Seq_Final),
    
    # --- Step 3: combined flag ---
    # If Ref or Alt is repetitive, mark as low-complexity
    Is_Low_Complexity = Ref_is_Repeat | Alt_is_Repeat
  ) %>% group_by(class) %>%
  summarise(
    total_count = n(),                          # Total sites in family
    repeat_count = sum(Is_Low_Complexity),      # Repetitive sites (TRUE count)
    repeat_rate = mean(Is_Low_Complexity)       # Repeat rate (TRUE proportion)
  ) %>%
  ungroup()

ggplot(plot_data, aes(x = reorder(class, repeat_rate), y = repeat_rate)) +
  
  # A. Bars
  geom_col(aes(fill = repeat_rate), width = 0.7) +
  
  # B. Text labels (percent + sample size)
  # label format: "75.0% (n=20)"
  geom_text(aes(label = paste0(sprintf("%.1f", repeat_rate * 100), "%", 
                               "\n(n=", total_count, ")")), 
            hjust = -0.1,  # Place text to the right of bars
            size = 3.5, 
            color = "black") +
  
  # C. Flip coordinates (long labels)
  coord_flip() +
  
  # D. Y-axis range (leave room for labels)
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0, max(plot_data$repeat_rate) * 1.2)) +
  
  # E. Color gradient (darker = higher repeat rate)
  scale_fill_gradient(low = "#85D4E3", high = "#F4B5BD") +
  
  # F. Labels and theme
  labs(
    title = "Proportion of Repetitive/Low-Complexity Sites by TF Family",
    x = "",  # Y axis already shows categories; X label can be empty
    y = "Percentage of Repetitive Sequences"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",              # Hide legend
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5) # Center title
  )
