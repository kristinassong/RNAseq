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


rule multiqc:
    input:
        fastqc_pretrim = expand(rules.fastqc_pretrim.output,sample=SAMPLES,read=['R1','R2']),
	fastqc_posttrim = expand(rules.fastqc_posttrim.output.html,sample=SAMPLES,read=['R1','R2']),
	picard = expand(rules.picard.output.pdf,sample=SAMPLES),
	kallisto = expand(rule.kallisto_quant.output,sample=SAMPLES)
    output:
        html = "results/multiqc/multiqc_report.html"
    params:
        scan_dir = "results/fastqc/pretrim results/fastqc/postrim results/STAR/*/*Log.final.out results/logs/kallisto results/picard results/logs/trimmomatic",
        outdir = directory("results/multiqc")
    log:
        "results/logs/multiqc.log"
    conda:
        "../envs/fastqc.yaml"
    message:
        "Summarize analysis results for multiple tools and samples in a single report."
    shell:
        "multiqc {params.scan_dir} "
        "--outdir {params.outdir} "
        "&> {log}"
