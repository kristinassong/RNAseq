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
        filtered_genes = rules.merge_kallisto_quant.output.tpm
    output:
        volcano = "results/volcano_plot/{comp}.svg"
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