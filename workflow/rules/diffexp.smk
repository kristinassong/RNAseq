### Adapted from Danny Bergeron's code

rule deseq2:
    input:
        quant = expand(rules.kallisto_quant.output, sample=SAMPLES),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = rules.tx2gene.output.tsv
    output:
        results = directory("results/deseq2"),
        out_files = comparisons_full
    params:
        kallisto_dir = "results/kallisto"
    log:
        "results/logs/deseq2.log"
    conda:
        "../envs/deseq2.yaml"
    message:
        "Perform differential expression analysis for various conditions."
    script:
        "../scripts/DESeq2_kallisto_tximport.R"

rule volcano_plot:
    input:
        DE_outdir = rules.deseq2.output.results,
        DE_output = "results/deseq2/{comp}.csv",
        filtered_genes = rules.merge_kallisto_quant.output.tpm,
        gtf = rules.filter_gtf_pc_genes.output.pc_gtf
    output:
        volcano = "results/volcano_plot/{comp}.svg",
        up_genes = "results/volcano_plot/{comp}_sig_diffexp_genes_up.tsv",
        down_genes = "results/volcano_plot/{comp}_sig_diffexp_genes_down.tsv"
    params:
        pval_threshold = 0.05
    log:
        "results/logs/volcano_plot/{comp}.log"
    conda:
        "../envs/python_plots.yaml"
    message:
        "Create a volcano plot using deseq2 output for each comparison in comparisons.tsv."
    script:
        "../scripts/volcano_plot.py"

rule GO_analysis_upregulated_genes:
    input:
        genes = rules.volcano_plot.output.up_genes
    output:
        bar_chart = "results/GO/GO_{comp}_upregulated_genes.svg"
    log:
        "results/logs/GO/{comp}_upregulated_genes.log"
    conda: 
        "../envs/GO.yaml"
    message:
        "GO analysis of upregulated genes in {wildcards.comp} represented as a bar chart."
    script:
        "../scripts/GO_bar_charts.py"

rule GO_analysis_downregulated_genes:
    input:
        genes = rules.volcano_plot.output.down_genes
    output:
        bar_chart = "results/GO/GO_{comp}_downregulated_genes.svg"
    log:
        "results/logs/GO/{comp}_downregulated_genes.log"
    conda: 
        "../envs/GO.yaml"
    message:
        "GO analysis of downregulated genes in {wildcards.comp} represented as a bar chart."
    script:
        "../scripts/GO_bar_charts.py"