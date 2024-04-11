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
    conda:
        "../envs/deseq2.yaml"
    message:
        "Perform differential expression analysis for {wildcards.comp}."
    script:
        "../scripts/DESeq2_kallisto_tximport.R"


rule volcano_plot:
    input:
        DE_outdir = rules.deseq2.output.results,
        DE_output = "results/deseq2/{comp}.csv",
        filtered_genes = rules.merge_kallisto_quant.output.tpm,
    output:
        volcano = "results/deseq2/{comp}.svg",
        up_genes = "results/deseq2/{comp}_sig_DE_up.tsv",
        down_genes = "results/deseq2/{comp}_sig_DE_down.tsv"
    params:
        pval_threshold = 0.05,
        gtf = config["path"]["genome_gtf"]
    conda:
        "../envs/GO.yaml"
    message:
        "Create a volcano plot using deseq2 output for {wildcards.comp}."
    script:
        "../scripts/volcano_plot.py"


rule GO_upregulated_genes:
    input:
        genes = rules.volcano_plot.output.up_genes
    output:
        bar_chart = "results/GO/GO_{comp}_up.svg"
    conda: 
        "../envs/GO.yaml"
    message:
        "GO analysis of upregulated genes in {wildcards.comp} represented as a bar chart."
    script:
        "../scripts/GO_bar_charts.py"


rule GO_downregulated_genes:
    input:
        genes = rules.volcano_plot.output.down_genes
    output:
        bar_chart = "results/GO/GO_{comp}_down.svg"
    conda: 
        "../envs/GO.yaml"
    message:
        "GO analysis of downregulated genes in {wildcards.comp} represented as a bar chart."
    script:
        "../scripts/GO_bar_charts.py"