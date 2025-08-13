
# DIFFERENTIAL ANALYSIS

# Raw counts and mapping statistics were used to calculate spike-in normalization factors

# Load necessary libraries 

library(edgeR); library(limma); library(dplyr)



# Load the spike-in normalized count table 

new_countdata_table <- read.table("raw_counts_si_normalized.txt", header = TRUE, sep = "\t")


# Set significance threshold for p-value

alpha <- 0.05


# Define experimental conditions for each sample in the count table

sample_info <- data.frame(
  row.names = colnames(new_countdata_table),
  condition = factor(rep(c("Input_KO_macroH2A1","Input_KO_rpb1","Input_KO_spt16",
                           "KO_macroH2A1","KO_rpb1","KO_spt16",
                           "Input_WT_macroH2A1","Input_WT_rpb1","Input_WT_spt16",
                           "WT_macroH2A1","WT_rpb1","WT_spt16"), each = 2)))


# Filter low read count genes

group <- factor(sample_info$condition)
dge <- DGEList(counts = new_countdata_table, group = sample_info$condition)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

write.table(dge$counts, "si_normalized_rawcounts_lowfiltered.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Differential analysis will compare knockout vs wildtype samples for:
# - macroH2A1 (KO_MACROH2A1 vs WT_MACROH2A1)
# - rpb1 (KO_RPB1 vs WT_RPB1)
# - spt16 (KO_SPT16 vs WT_SPT16)
#
# Additionally, comparisons for Input samples will be performed:
# - macroH2A1 (INPUT_KO_MACROH2A1 vs INPUT_WT_MACROH2A1)
# - rpb1 (INPUT_KO_RPB1 vs INPUT_WT_RPB1)
# - spt16 (INPUT_KO_SPT16 vs INPUT_WT_SPT16)


# ------- KO_MACROH2A1 VS WT_MACROH2A1 ------- 


# Design condition matrix

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition) 
contrast_matrix <- makeContrasts(KOvsWT_macroH2A1 = KO_macroH2A1 - WT_macroH2A1, levels = design)


# Differential analysis using limma with voom transformation

dge_voom <- voom(dge, design, plot = FALSE)
fit_limma <- lmFit(dge_voom, design)
fit_limma <- contrasts.fit(fit_limma, contrast_matrix)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = "KOvsWT_macroH2A1", number = nrow(dge_voom$E), adjust.method = "BH")
sig_limma <- sum(res_limma$P.Value < alpha)
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_macroH2A1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_macroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_macroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)



# ------- INPUT_KO_MACROH2A1 vs INPUT_WT_MACROH2A1 ------- 

# Design condition matrix

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition) 
contrast_matrix <- makeContrasts(KOvsWT_InputmacroH2A1 = Input_KO_macroH2A1 - Input_WT_macroH2A1, levels = design)


# Differential analysis using limma with voom transformation

dge_voom <- voom(dge, design, plot = FALSE)
fit_limma <- lmFit(dge_voom, design)
fit_limma <- contrasts.fit(fit_limma, contrast_matrix)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = "KOvsWT_InputmacroH2A1", number = nrow(dge_voom$E), adjust.method = "BH")
sig_limma <- sum(res_limma$P.Value < alpha)
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_InputmacroH2A1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_InputmacroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_InputmacroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ------- KO_RPB1 VS WT_RPB1 ------- 


# Design condition matrix

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition) 
contrast_matrix <- makeContrasts(KOvsWT_rpb1 = KO_rpb1 - WT_rpb1, levels = design)


# Differential analysis using limma with voom transformation

dge_voom <- voom(dge, design, plot = FALSE)
fit_limma <- lmFit(dge_voom, design)
fit_limma <- contrasts.fit(fit_limma, contrast_matrix)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = "KOvsWT_rpb1", number = nrow(dge_voom$E), adjust.method = "BH")
sig_limma <- sum(res_limma$P.Value < alpha)
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_rpb1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_rpb1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_rpb1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ------- INPUT_KO_RPB1 vs INPUT_WT_RPB1 ------- 


# Design condition matrix
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition) 
contrast_matrix <- makeContrasts(KOvsWT_Inputrpb1 = Input_KO_rpb1 - Input_WT_rpb1, levels = design)


# Differential analysis using limma with voom transformation
dge_voom <- voom(dge, design, plot = FALSE)
fit_limma <- lmFit(dge_voom, design)
fit_limma <- contrasts.fit(fit_limma, contrast_matrix)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = "KOvsWT_Inputrpb1", number = nrow(dge_voom$E), adjust.method = "BH")
sig_limma <- sum(res_limma$P.Value < alpha)
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching
expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_Inputrpb1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_Inputrpb1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_Inputrpb1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ------- KO_SPT16 VS WT_SPT16 ------- 


# Design condition matrix
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition) 
contrast_matrix <- makeContrasts(KOvsWT_spt16 = KO_spt16 - WT_spt16, levels = design)


# Differential analysis using limma with voom transformation
dge_voom <- voom(dge, design, plot = FALSE)
fit_limma <- lmFit(dge_voom, design)
fit_limma <- contrasts.fit(fit_limma, contrast_matrix)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = "KOvsWT_spt16", number = nrow(dge_voom$E), adjust.method = "BH")
sig_limma <- sum(res_limma$P.Value < alpha)
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_spt16_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL 


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_spt16.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_spt16.txt", sep = "\t", quote = FALSE, row.names = FALSE)



# ------- INPUT_KO_SPT16 vs INPUT_WT_SPT16 ------- 


# Design condition matrix

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition) 
contrast_matrix <- makeContrasts(KOvsWT_Inputspt16 = Input_KO_spt16 - Input_WT_spt16, levels = design)


# Differential analysis using limma with voom transformation

dge_voom <- voom(dge, design, plot = FALSE)
fit_limma <- lmFit(dge_voom, design)
fit_limma <- contrasts.fit(fit_limma, contrast_matrix)
fit_limma <- eBayes(fit_limma)
res_limma <- topTable(fit_limma, coef = "KOvsWT_Inputspt16", number = nrow(dge_voom$E), adjust.method = "BH")
sig_limma <- sum(res_limma$P.Value < alpha)
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_Inputspt16_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_Inputspt16.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_Inputspt16.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# -----------------------------------------------------------
# ChIP-seq sample ratio calculation and differential analysis
  
  
# Load normalized and filtered counts, then compute ratios for replicates 1 and 2.
  
# ------- SPT16/MACROH2A1 RATIO ------- 

normcounts <- read.delim("si_normalized_rawcounts_lowfiltered.txt", sep = "\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)

ratiospt16macroH2A1 <- data.frame(
  KO_rep1 = normcounts$PF5.KO_spt16_rep1.bam / normcounts$PF5.KO_macroH2A1_rep1.bam,
  KO_rep2 = normcounts$PF5.KO_spt16_rep2.bam / normcounts$PF5.KO_macroH2A1_rep2.bam,
  WT_rep1 = normcounts$WT_spt16_rep1.bam / normcounts$WT_macroH2A1_rep1.bam,
  WT_rep2 = normcounts$WT_spt16_rep2.bam / normcounts$WT_macroH2A1_rep2.bam,
  row.names = rownames(normcounts)
)

ratiospt16macroH2A1 <- as.data.frame(ratiospt16macroH2A1)
ratiospt16macroH2A1[is.na(ratiospt16macroH2A1)] <- 0
ratiospt16macroH2A1[!is.finite(as.matrix(ratiospt16macroH2A1))] <- 0
sum(is.na(ratiospt16macroH2A1))                       
sum(is.infinite(as.matrix(ratiospt16macroH2A1)))      
write.table (ratiospt16macroH2A1, "ratio_spt16macroH2A1.txt")


# Differential analysis 

alpha <- 0.05 


# Preprocessing: log2 transformation with pseudocount and cleanup of non-finite values

log2_ratios <- log2(ratiospt16macroH2A1 + 1e-6)  
log2_ratios[!is.finite(as.matrix(log2_ratios))] <- 0  


# Design condition matrix

sample_info <- data.frame(
  row.names = colnames(log2_ratios),
  condition = factor(rep(c("KO", "WT"), each = 2))
)
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)
contrast_matrix <- makeContrasts(KOvsWT_ratio = KO - WT, levels = design)


# Limma differential analysis

fit <- lmFit(log2_ratios, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef = "KOvsWT_ratio", number = nrow(log2_ratios), adjust.method = "BH")
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]
cat("Differential genes (p < ", alpha, "): ", nrow(res_limma_sig), "\n")


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_ratio_spt16macroH2A1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_ratio_spt16macroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_ratio_spt16macroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ------- SPT16/RPB1 RATIO ------- 

normcounts <- read.delim("si_normalized_rawcounts_lowfiltered.txt", sep = "\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)

ratiospt16rpb1 <- data.frame(
  KO_rep1 = normcounts$PF5.KO_spt16_rep1.bam / normcounts$PF5.KO_Rpb1_rep1.bam,
  KO_rep2 = normcounts$PF5.KO_spt16_rep2.bam / normcounts$PF5.KO_Rpb1_rep2.bam,
  WT_rep1 = normcounts$WT_spt16_rep1.bam / normcounts$WT_Rpb1_rep1.bam,
  WT_rep2 = normcounts$WT_spt16_rep2.bam / normcounts$WT_Rpb1_rep2.bam,
  row.names = rownames(normcounts)
)

ratiospt16rpb1 <- as.data.frame(ratiospt16rpb1)
ratiospt16rpb1[is.na(ratiospt16rpb1)] <- 0
ratiospt16rpb1[!is.finite(as.matrix(ratiospt16rpb1))] <- 0
sum(is.na(ratiospt16rpb1))                        
sum(is.infinite(as.matrix(ratiospt16rpb1)))       
write.table (ratiospt16rpb1, "ratiospt16rpb1.txt")


# Differential analysis 

alpha <- 0.05  


# Preprocessing: log2 transformation with pseudocount and cleanup of non-finite values

log2_ratios <- log2(ratiospt16rpb1 + 1e-6)  
log2_ratios[!is.finite(as.matrix(log2_ratios))] <- 0 


# Design condition matrix

sample_info <- data.frame(
  row.names = colnames(log2_ratios),
  condition = factor(rep(c("KO", "WT"), each = 2))
)
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)
contrast_matrix <- makeContrasts(KOvsWT_ratio = KO - WT, levels = design)


# Limma differential analysis

fit <- lmFit(log2_ratios, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef = "KOvsWT_ratio", number = nrow(log2_ratios), adjust.method = "BH")
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]
cat("Differential genes (p < ", alpha, "): ", nrow(res_limma_sig), "\n")


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_ratio_spt16rpb1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_ratio_spt16rpb1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_ratio_spt16rpb1.txt", sep = "\t", quote = FALSE, row.names = FALSE)




 # ------- RPB1/MACROH2A1 RATIO ------- 


normcounts <- read.delim("si_normalized_rawcounts_lowfiltered.txt", sep = "\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
ratiorpb1macroh2a1 <- data.frame(
  KO_rep1 = normcounts$PF5.KO_Rpb1_rep1.bam / normcounts$PF5.KO_macroH2A1_rep1.bam,
  KO_rep2 = normcounts$PF5.KO_Rpb1_rep2.bam / normcounts$PF5.KO_macroH2A1_rep2.bam,
  WT_rep1 = normcounts$WT_Rpb1_rep1.bam / normcounts$WT_macroH2A1_rep1.bam,
  WT_rep2 = normcounts$WT_Rpb1_rep2.bam / normcounts$WT_macroH2A1_rep2.bam,
  row.names = rownames(normcounts)
)
ratiorpb1macroh2a1 <- as.data.frame(ratiorpb1macroh2a1)
ratiorpb1macroh2a1[is.na(ratiorpb1macroh2a1)] <- 0
ratiorpb1macroh2a1[!is.finite(as.matrix(ratiorpb1macroh2a1))] <- 0
sum(is.na(ratiorpb1macroh2a1))                        
sum(is.infinite(as.matrix(ratiorpb1macroh2a1)))       
write.table (ratiorpb1macroh2a1, "ratiorpb1macroH2A1.txt")


# Differential analysis 

alpha <- 0.05  


# Preprocessing: log2 transformation with pseudocount and cleanup of non-finite values

log2_ratios <- log2(ratiorpb1macroh2a1 + 1e-6)  
log2_ratios[!is.finite(as.matrix(log2_ratios))] <- 0  


# Design condition matrix

sample_info <- data.frame(
  row.names = colnames(log2_ratios),
  condition = factor(rep(c("KO", "WT"), each = 2))
)
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)
contrast_matrix <- makeContrasts(KOvsWT_ratio = KO - WT, levels = design)


# Limma differential analysis

fit <- lmFit(log2_ratios, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef = "KOvsWT_ratio", number = nrow(log2_ratios), adjust.method = "BH")
res_limma_sig <- res_limma[res_limma$P.Value < alpha, ]
cat("Differential genes (p < ", alpha, "): ", nrow(res_limma_sig), "\n")


# Intersect with expressed genes from original .bed file. 'GENEID' column was 
# removed and set as row names for matching

expressed_genes <- read.table("NOGENEID_HCT116_expressed_genes_SI_computematrix.txt", row.names= 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_limma_sig_filtered <- res_limma_sig[rownames(res_limma_sig) %in% rownames(expressed_genes), ]
write.table(res_limma_sig_filtered, "limma_results_ratio_rpb1macroH2A1_KOvsWT_filtered.txt", sep = "\t", quote = FALSE)
res_limma_sig_filtered$GeneID <- rownames(res_limma_sig_filtered)
rownames(res_limma_sig_filtered) <- NULL  


# Filter and save genes with lower signal in KO (logFC < -0.5)

res_limma_sig_filtered_down <- res_limma_sig_filtered %>%
  filter(logFC < -0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_down, "DOWNDGS_limma_ratio_rpb1macroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Filter and save genes with higher signal in KO (logFC > 0.5)

res_limma_sig_filtered_up <- res_limma_sig_filtered %>%
  filter(logFC > 0.5) %>%
  select(GeneID, logFC)
write.table(res_limma_sig_filtered_up, "UPDGS_limma_ratio_rpb1macroH2A1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


