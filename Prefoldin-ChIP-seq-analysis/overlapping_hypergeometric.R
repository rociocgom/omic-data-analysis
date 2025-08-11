
# OVERLAPPINGS AND HYPERGEOMETRIC TEST


# Load libraries

library(tidyverse)


# Load all files: define a function to read only the first column (genes) 

read_genes <- function(path) {
  read_tsv(path, skip = 1, col_names = FALSE, show_col_types = FALSE)[[1]]
}

# Load gene lists by extracting only the first column

gene_lists <- list(
  #rpb1
  rpb1_DOWN = read_genes("DEF/DIF_ANALYSIS/RPB1_RESULTS/DOWNDEGS_limma_rpb1.txt"),
  rpb1_UP = read_genes("DEF/DIF_ANALYSIS/RPB1_RESULTS/UPDEGS_limma_rpb1.txt"),
  
  #spt16
  spt16_DOWN = read_genes("DEF/DIF_ANALYSIS/SPT16_RESULTS/DOWNDEGS_limma_spt16.txt"),
  spt16_UP = read_genes("DEF/DIF_ANALYSIS/SPT16_RESULTS/UPDEGS_limma_spt16.txt"),
  
  #macroh2a1
  macroH2A1_DOWN = read_genes("DEF/DIF_ANALYSIS/MACROH2A1_RESULTS/DOWNDEGS_limma_macroH2A1.txt"),
  macroH2A1_UP = read_genes("DEF/DIF_ANALYSIS/MACROH2A1_RESULTS/UPDEGS_limma_macroH2A1.txt"),
  
  #ratio spt16-macroH2A1
  ratio_spt16_macro_DOWN = read_genes("DEF/DIF_ANALYSIS/RATIO SPT16_MACROH2A1 RESULTS/DOWNDEGS_limma_ratio_spt16macroH2A1.txt"),
  ratio_spt16_macro_UP = read_genes("DEF/DIF_ANALYSIS/RATIO SPT16_MACROH2A1 RESULTS/UPDEGS_limma_ratio_spt16macroH2A1.txt"),
  
  #ratio rpb1-macroH2A1
  ratio_rpb1_macro_DOWN = read_genes("DEF/DIF_ANALYSIS/RATIO RPB1_MACROH2A1 RESULTS/DOWNDEGS_limma_ratio_rpb1macroH2A1.txt"),
  ratio_rpb1_macro_UP = read_genes("DEF/DIF_ANALYSIS/RATIO RPB1_MACROH2A1 RESULTS/UPDEGS_limma_ratio_rpb1macroH2A1.txt"),
  
  #ratio spt16-rpb1
  ratio_spt16_rpb1_DOWN = read_genes("DEF/DIF_ANALYSIS/RATIO SPT16_RPB1 RESULTS/DOWNDEGS_limma_ratio_spt16rpb1.txt"),
  ratio_spt16_rpb1_UP = read_genes("DEF/DIF_ANALYSIS/RATIO SPT16_RPB1 RESULTS/UPDEGS_limma_ratio_spt16rpb1.txt")
)


# Associate each list with its base condition

conds <- names(gene_lists)

elements <- c(
  rpb1_DOWN             = "rpb1",
  rpb1_UP               = "rpb1",
  spt16_DOWN            = "spt16",
  spt16_UP              = "spt16",
  macroH2A1_DOWN        = "macroH2A1",
  macroH2A1_UP          = "macroH2A1",
  
  ratio_spt16_macro_DOWN = "spt16_macro",
  ratio_spt16_macro_UP   = "spt16_macro",
  ratio_spt16_rpb1_DOWN  = "spt16_rpb1",
  ratio_spt16_rpb1_UP    = "spt16_rpb1",
  ratio_rpb1_macro_DOWN  = "rpb1_macro",
  ratio_rpb1_macro_UP    = "rpb1_macro"
)


# Create non-redundant pairwise comparisons across different conditions

pairs <- expand.grid(a = conds, b = conds, stringsAsFactors = FALSE) %>%
  filter(a < b) %>%      
  filter(elements[a] != elements[b])   
nrow(pairs)


# Hypergeometric overlap analysis for all gene list pairs

results <- tibble()  

N_total <- 12763  

for(i in seq_len(nrow(pairs))) {
  la <- pairs$a[i]
  lb <- pairs$b[i]
  va <- gene_lists[[la]]
  vb <- gene_lists[[lb]]
  m <- length(va)      
  n <- N_total - m    
  k <- length(vb)          
  genes_overlap <- intersect(va, vb)
  x <- length(intersect(va, vb)) 
  pval <- phyper(q = x - 1, m = m, n = n, k = k, lower.tail = FALSE)

  genes_str <- if (x > 0) paste(genes_overlap, collapse = ";") else NA_character_
  
# Save results for each comparison
  results <- bind_rows(results,
                       tibble(
                         list1    = la,
                         list2    = lb,
                         n1       = m,
                         n2       = k,
                         overlap  = x,
                         not_in_list1 = m - x,  # Genes in A not overlapping
                         not_in_list2 = k - x,  # # Genes in B not overlapping

                         pvalue   = pval,
                         genes = genes_str
                       ))
}



# Adjust p-values for multiple testing, filter significant overlaps, and save results

results <- results %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(padj)

sig_pairs <- results %>% filter(padj < 0.05)

write.table (results, "overlap_results_with_genes.txt", sep = "\t", row.names = FALSE)
write.table (sig_pairs, "differential_overlap_results_with_genes.txt", sep = "\t", row.names = FALSE)

