
# ==============================================================================
# Bulk RNA-seq Analysis with DESeq2
# ==============================================================================

# Packages ----------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Set working directory
setwd("/Users/hyjung/Gladstone Dropbox/Hyuncheol Jung/20260128BulkRNAseq3rd")

# Create output directories
dir.create("output", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# Import data -------------------------------------------------------------
counts <- read_csv("combined_counts_nostim.csv", show_col_types = FALSE)
counts

# Check column names
colnames(counts)
# [1] "gene_symbol" "ensembl_id"  "ntd_d1" "aavs_d1" "a20_d1" "ntd_d2" ...

# Prepare data for DESeq2 -------------------------------------------------

# Count matrix format for DESeq2 input
# Remove gene_symbol, use ensembl_id as rownames
cts <- counts %>%
  select(-gene_symbol) %>%
  distinct(ensembl_id, .keep_all = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("ensembl_id") %>%
  as.matrix()

# DESeq2 requires integer counts - round if needed
cts <- round(cts)

head(cts)
dim(cts)  # Check dimensions

## Save DESeq2 input
cts %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  write_tsv("output/DESeq2_input_cts.txt")

# Build condition table ---------------------------------------------------
# Your samples: ntd_d1, aavs_d1, a20_d1, ntd_d2, aavs_d2, a20_d2, ... a20_d6
# Format: {sgrna}_{donor}

condition_df <- tibble(Sample = colnames(cts)) %>%
  mutate(
    # Extract sgRNA type (ntd, aavs, a20)
    sgrna = str_extract(Sample, "^[^_]+"),
    sgrna = case_when(
      sgrna == "ntd" ~ "NTD",
      sgrna == "aavs" ~ "AAVS1", 
      sgrna == "a20" ~ "A20",
      TRUE ~ sgrna
    ),
    # Extract donor (d1, d2, d3, d4, d5, d6)
    donor = str_extract(Sample, "d[0-9]+"),
    donor = toupper(donor)
  )

condition_df
# Convert to factors
condition_df <- condition_df %>%
  mutate(
    sgrna = factor(sgrna, levels = c("NTD", "AAVS1", "A20")),
    donor = factor(donor)
  )

condition_df

# Create coldata for DESeq2
coldata <- condition_df %>%
  as.data.frame() %>%
  column_to_rownames("Sample")

coldata

# Check that the two dataframes are in agreement
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

## Save coldata
coldata %>%
  rownames_to_column("Sample") %>%
  write_tsv("output/DESeq2_coldata.txt")

# DESeq2 pipeline ---------------------------------------------------------

# Prepare DESeq2 - design with donor as batch effect
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ donor + sgrna  # donor as batch, sgrna as main effect
)

dds

# Filter low count genes (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Run DESeq
dds <- DESeq(dds)
dds

# Check results names
resultsNames(dds)

# Get variance stabilized transformed data (normalized data)
vsd <- vst(dds, blind = FALSE)

## Save vsd
assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  write_tsv("output/DESeq2_varianceStabilizedTransformedData.txt")

# Global analysis ---------------------------------------------------------

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

cor_heatmap <- pheatmap(
  cor(assay(vsd)),
  col = viridisLite::viridis(100),
  main = "Sample Pearson Correlations",
  annotation_col = coldata[, c("sgrna", "donor"), drop = FALSE]
)

ggsave("plots/01.01_correlation-heatmap.png", 
       plot = cor_heatmap,
       width = 8, height = 7)

# PCA ---------------------------------------------------------------------
pcaData <- plotPCA(vsd, intgroup = c("sgrna", "donor"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- pcaData %>%
  ggplot(aes(x = PC1, y = PC2, color = sgrna, shape = donor, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot - Bulk RNA-seq")

pca_plot

ggsave("plots/01.02_PCA-plot.png", 
       plot = pca_plot,
       width = 8, height = 6)

# Differential Expression Analysis ----------------------------------------

# Compare A20 vs AAVS1 (control)
res_A20_vs_AAVS1 <- results(dds, contrast = c("sgrna", "A20", "AAVS1"))
res_A20_vs_AAVS1 <- res_A20_vs_AAVS1[order(res_A20_vs_AAVS1$pvalue), ]
summary(res_A20_vs_AAVS1)

# Save results
res_A20_vs_AAVS1 %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(counts %>% select(ensembl_id, gene_symbol), by = "ensembl_id") %>%
  relocate(gene_symbol, .after = ensembl_id) %>%
  arrange(pvalue) %>%
  write_csv("output/DESeq2_results_A20_vs_AAVS1.csv")

# Compare A20 vs NTD
res_A20_vs_NTD <- results(dds, contrast = c("sgrna", "A20", "NTD"))
res_A20_vs_NTD <- res_A20_vs_NTD[order(res_A20_vs_NTD$pvalue), ]
summary(res_A20_vs_NTD)

# Save results
res_A20_vs_NTD %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(counts %>% select(ensembl_id, gene_symbol), by = "ensembl_id") %>%
  relocate(gene_symbol, .after = ensembl_id) %>%
  arrange(pvalue) %>%
  write_csv("output/DESeq2_results_A20_vs_NTD.csv")

# Volcano Plot ------------------------------------------------------------
# For A20 vs AAVS1

volcano_df <- res_A20_vs_AAVS1 %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(counts %>% select(ensembl_id, gene_symbol), by = "ensembl_id") %>%
  mutate(
    group = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "Up",
      log2FoldChange < -1 & padj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano_plot <- volcano_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#D7191C", "Down" = "#2C7BB6", "NS" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = filter(volcano_df, group %in% c("Up", "Down")) %>% slice_min(padj, n = 20),
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    title = "Volcano Plot: A20 vs AAVS1",
    x = expression(log[2]~FoldChange),
    y = expression(-log[10]~padj)
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

volcano_plot

ggsave("plots/01.03_volcano-A20-vs-AAVS1.png",
       plot = volcano_plot,
       width = 8, height = 6)



# PCA - using PC2 and PC3 -----------------------------------------------------

# Get all PCs (not just PC1 and PC2)
pcaData_full <- prcomp(t(assay(vsd)))

# Calculate percent variance for each PC
percentVar <- round(100 * (pcaData_full$sdev^2 / sum(pcaData_full$sdev^2)))

# Create dataframe with PC2 and PC3
pcaData <- data.frame(
  PC2 = pcaData_full$x[, 2],
  PC3 = pcaData_full$x[, 3],
  name = rownames(pcaData_full$x)
) %>%
  left_join(condition_df, by = c("name" = "Sample"))

pca_plot <- pcaData %>%
  ggplot(aes(x = PC2, y = PC3, color = sgrna, shape = donor, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  xlab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylab(paste0("PC3: ", percentVar[3], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot (PC2 vs PC3)")

pca_plot

ggsave("plots/01.02_PCA-plot_PC2-PC3.png", 
       plot = pca_plot,
       width = 8, height = 6)






# SessionInfo -------------------------------------------------------------
sessionInfo()
