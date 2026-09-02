import os

# Automatically detect samples from the data/raw/ directory
SAMPLES, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")

# Target output: We want the final DE results and the QC reports
rule all:
    input:
        expand("results/qc/{sample}_fastp.html", sample=SAMPLES),
        "results/DE_results/significant_genes.csv",
        "results/DE_results/volcano_plot_summary.png"

rule fastp_trimming:
    input:
        r1 = "data/raw/{sample}_R1.fastq.gz",
        r2 = "data/raw/{sample}_R2.fastq.gz"
    output:
        r1 = "data/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2 = "data/trimmed/{sample}_trimmed_R2.fastq.gz",
        html = "results/qc/{sample}_fastp.html"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --cut_front --cut_tail --thread {threads} --html {output.html}"

rule salmon_quant:
    input:
        r1 = "data/trimmed/{sample}_trimmed_R1.fastq.gz",
        r2 = "data/trimmed/{sample}_trimmed_R2.fastq.gz",
        index = "references/salmon_index"
    output:
        quant_dir = directory("quants/{sample}_quant"),
        sf = "quants/{sample}_quant/quant.sf"
    threads: 8
    conda:
        "envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input.index} \
                     -l A \
                     -1 {input.r1} -2 {input.r2} \
                     --validateMappings \
                     --gcBias \
                     --seqBias \
                     --minAssignedFrags 1 \
                     -p {threads} \
                     -o {output.quant_dir}
        """

rule deseq2_analysis:
    input:
        quants = expand("quants/{sample}_quant/quant.sf", sample=SAMPLES),
        tx2gene = "references/tx2gene.csv"
    output:
        results_csv = "results/DE_results/significant_genes.csv"
    conda:
        "envs/deseq2.yaml"
    shell:
        "Rscript scripts/run_DE.R {input.tx2gene} {output.results_csv}"

rule volcano_plot:
    input:
        results_csv = "results/DE_results/significant_genes.csv"
    output:
        plot = "results/DE_results/volcano_plot_summary.png"
    conda:
        "envs/deseq2.yaml"
    shell:
        "Rscript scripts/volcano_plot.R"