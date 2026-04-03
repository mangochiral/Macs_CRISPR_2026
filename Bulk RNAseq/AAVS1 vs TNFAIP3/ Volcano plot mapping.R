# Get genes from specific Hallmark pathways
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# Define pathways to highlight (use exact names from MSigDB)
pathways_to_highlight <- list(
  "Inflammatory Response" = hallmark_sets %>% 
    filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% 
    pull(gene_symbol),
  
  "TNFa Signaling via NFkB" = hallmark_sets %>% 
    filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>% 
    pull(gene_symbol)
)

# Check how many genes in each pathway
lapply(pathways_to_highlight, length)

# Create volcano plot dataframe with pathway annotations
volcano_df <- de_results %>%
  mutate(
    group = case_when(
      log2FoldChange > 2 & padj < 0.05 ~ "Up",
      log2FoldChange < -2 & padj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Add pathway membership
pathway_colors <- c("#E41A1C", "#377EB8")
pathway_names <- names(pathways_to_highlight)

volcano_df$pathway <- "None"
for (i in seq_along(pathways_to_highlight)) {
  pathway_genes <- pathways_to_highlight[[i]]
  volcano_df$pathway[volcano_df$gene_symbol %in% pathway_genes] <- pathway_names[i]
}

# Create factor with "None" last for plotting order
volcano_df$pathway <- factor(volcano_df$pathway, 
                             levels = c(pathway_names, "None"))

# Volcano plot with pathway highlighting
# Volcano plot - only show significant pathway genes
volcano_pathway_plot <- volcano_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  # Background points (non-pathway genes)
  geom_point(data = filter(volcano_df, pathway == "None"),
             color = "grey80", size = 0.1, alpha = 0.5) +
  # Pathway genes - ONLY SIGNIFICANT
  geom_point(data = filter(volcano_df, pathway != "None" & padj < 0.05 & abs(log2FoldChange) > 1),
             aes(color = pathway), size = 1, alpha = 0.8) +
  # Color scale for pathways
  scale_color_manual(values = setNames(pathway_colors, pathway_names),
                     name = "Pathway") +
  # Reference lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Label significant pathway genes
  geom_text_repel(
    data = filter(volcano_df, pathway != "None" & padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene_symbol, color = pathway),
    size = 2,
    box.padding = 0.3,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  # Axis limits
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 10)) +
  labs(
    title = "Volcano Plot: A20 vs AAVS1",
    subtitle = "Significant Hallmark pathway genes highlighted",
    x = expression(log[2]~FoldChange),
    y = expression(-log[10]~padj)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

volcano_pathway_plot

ggsave("plots/02.04_volcano_with_hallmark_pathways.png",
       plot = volcano_pathway_plot,
       width = 10, height = 8)
