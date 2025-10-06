# Differential Gene Expression using Deseq2

# Load required libraries
library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(apeglm)
library(org.Mm.eg.db)  # Mouse genome database

# Read the gene expression matrix
counts_matrix <- read.table("gene_expression_matrix.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Clean up column names - remove path and keep only sample names
colnames(counts_matrix) <- gsub("RSEM_results/", "", colnames(counts_matrix))
colnames(counts_matrix) <- gsub("\\.genes\\.results", "", colnames(counts_matrix))



# Create sample information
samdf <- data.frame(
  sample = colnames(counts_matrix),
  group = ifelse(grepl("DMSOR", colnames(counts_matrix)), "DMSO", "DOX")
)

# Ensure group is a factor with appropriate levels
samdf$group <- factor(samdf$group, levels = c("DMSO", "DOX"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),  # DESeq2 expects integer counts
  colData = samdf,
  design = ~ group
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("group", "DOX", "DMSO"))
res <- lfcShrink(dds, coef = "group_DOX_vs_DMSO", type = "apeglm")

# Convert to data frame and remove NA values
res_df <- as.data.frame(res) %>%
  na.omit() %>%
  mutate(
    gene = rownames(.),
    significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      padj < 0.05 & abs(log2FoldChange) <= 1 ~ "P-value only",
      padj >= 0.05 & abs(log2FoldChange) > 1 ~ "Fold change only",
      TRUE ~ "Not significant"
    )
  )

# Print summary
print("Differential expression summary:")
print(table(res_df$significance))

# Identify top significant genes 
top_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(20)

# Map ENSEMBL to gene symbols
ensembl_ids <- gsub("\\.[0-9]+$", "", rownames(res_df))

gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add gene symbols to top_genes
res_df$gene_symbol <- gene_symbols


# Create a new column for direction of significant genes
res_df <- res_df %>%
  mutate(
    direction = case_when(
      padj >= 0.05 | abs(log2FoldChange) <= 1 ~ "Not significant",
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated"
    )
  )
# Also add direction to top_genes for labeling
top_genes <- top_genes %>%
  mutate(
    direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")
  )

# Create volcano plot with red/blue colors
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Downregulated" = "blue",
    "Not significant" = "gray", 
    "Upregulated" = "red"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_symbol, color = direction),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.2,
    show.legend = FALSE  # Hide legend for labels
  ) +
  labs(
    title = "Volcano Plot: DOX vs DMSO (DESeq2)",
    x = "Log2 Fold Change (DOX/DMSO)",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom"
  )

# Print the plot
print(volcano_plot)
