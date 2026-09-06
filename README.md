cat << 'EOF' > README.md
# Automated End-to-End RNA-Seq Pipeline

A fully reproducible, automated bulk RNA-Seq analysis pipeline built using **Snakemake**, **Conda/Mamba**, and **R (DESeq2)**. This pipeline handles raw FASTQ inputs all the way to differential gene expression analysis and data visualization.

---

## 🚀 Pipeline Architecture
1. **Quality Control & Trimming**: Uses `fastp` to filter reads and generate interactive HTML QC reports.
2. **Transcript Quantification**: Uses `Salmon` for fast, lightweight pseudo-alignment against a GENCODE reference transcriptome index.
3. **Differential Expression Analysis**: Uses Bioconductor's `tximport` and `DESeq2` to process count matrices, handle transcript-level uncertainty, and calculate adjusted p-values ($padj < 0.05$).
4. **Data Visualization**: Automatically generates a publication-ready **Volcano Plot** summarizing up- and down-regulated genes.

---

## 🛠️ Tech Stack & Tools
* **Workflow Manager**: Snakemake (v9+)
* **Environment Management**: Conda / Mamba (Isolated YAML-based environments)
* **Alignment/Quantification**: Salmon
* * **QC**: fastp
* **Statistical Analysis**: R, DESeq2, tximport, jsonlite
* **Visualization**: base R graphics / ggplot-ready outputs

---

## 💻 Challenges Overcome During Development
* **Resource Management**: Handled Linux OOM (Out-Of-Memory) limits on an 8GB RAM laptop by optimizing core throttling and memory allocation during Salmon index loading.
* **Conda Environment Isolation**: Resolved package dependencies (adding `r-jsonlite`) to support Salmon inferential replicate imports in R.
* **ID Version Harmonization**: Handled GENCODE transcript-to-gene mapping discrepancies by writing custom regex handling inside R (`tximport` with custom suffix trimming).

---

## 📂 Project Structure
```text
RNA-Seq/
├── Snakefile                # Main Snakemake workflow orchestration
├── envs/                    # Isolated Conda YAML environments (salmon, deseq2, etc.)
├── scripts/                 # Custom R scripts for DESeq2 and Volcano plots
├── references/              # Genome indices and tx2gene.csv mapping files
├── data/raw/                # Input FASTQ files (Control & Treated)
└── results/                 # Final output metrics, counts, and volcano plots
