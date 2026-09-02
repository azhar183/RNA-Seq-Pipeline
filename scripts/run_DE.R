suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
tx2gene_path <- args[1]
output_path <- args[2]

quant_files <- list.files("quants", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
if(length(quant_files) == 0) {
    stop("No quant.sf files found in quants/ directory. Make sure you have fastq files in data/raw/")
}
sample_names <- gsub("_quant", "", basename(dirname(quant_files)))
names(quant_files) <- sample_names

conditions <- ifelse(grepl("Treated", sample_names, ignore.case=TRUE), "Treated", "Control")
colData <- data.frame(condition = factor(conditions))
rownames(colData) <- sample_names

tx2gene <- read_csv(tx2gene_path, show_col_types = FALSE)

# THE REAL FIX: Chop the decimals off the tx2gene transcript IDs so they match Salmon
tx2gene[[1]] <- gsub("\\..*", "", tx2gene[[1]])

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

sig_res <- subset(res, padj < 0.05)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(sig_res), file = output_path)
