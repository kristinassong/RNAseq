rule pca:
    input:
        tpm = rules.merge_kallisto_quant.output.tpm
    output:
        plot = "results/pca/pca_{comp}.svg",
        tsv = "results/pca/pca_{comp}.tsv"
    params:
        design = 'resources/design.tsv'
    conda:
        "../envs/plots.yaml"
    message:
        "Generate a PCA plot to observe variance between {wildcards.comp} samples."
    script:
        "../scripts/pca.py"