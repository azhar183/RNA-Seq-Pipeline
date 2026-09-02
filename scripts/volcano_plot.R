suppressPackageStartupMessages({
  library(readr)
})

# Read the significant genes CSV file
res_path <- "results/DE_results/significant_genes.csv"
if(!file.exists(res_path)) {
    stop("Results file not found. Make sure the pipeline has run successfully.")
}

res <- read_csv(res_path, show_col_types = FALSE)

# Basic console summary
cat("Successfully loaded results with", nrow(res), "significant genes.\n")
cat("Top regulated gene:\n")
print(res[1, c(1, 2, 3, 7)])

# If ggplot2 is available, generate a simple text/ascii summary or save a plot
png("results/DE_results/volcano_plot_summary.png", width=800, height=600)
plot(res$log2FoldChange, -log10(res$padj), 
     col = ifelse(res$log2FoldChange > 0, "red", "blue"),
     pch = 20, main = "RNA-Seq Volcano Plot (Test Data)",
     xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted P-value")
abline(h = -log10(0.05), col = "gray", lty = 2)
dev.off()

cat("Volcano plot saved to results/DE_results/volcano_plot_summary.png!\n")
