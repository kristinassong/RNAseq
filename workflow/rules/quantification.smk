rule build_transcriptome:
    input:
        genome = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        config["path"]["transcriptome"]
    conda:
        "../envs/gffread.yaml"
    message:
        "Build a reference transcriptome using gffread."
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output}"

rule kallisto_index:
    input:
        rules.build_transcriptome.output
    output:
        "data/references/kallisto.idx"
    params:
        31
    conda:
        "../envs/kallisto.yaml"
    log:
        "results/logs/kallisto/index.log"
    message:
        "Builds an index from the FASTA file."
    shell:
        "kallisto index "
        "--index={output} "
        "--kmer-size={params} "
        "{input} "
        "&> {log}"

rule kallisto_quant:
    input:
        idx = rules.kallisto_index.output,
        fq1 = rules.trimmomatic.output.r1,
        fq2 = rules.trimmomatic.output.r2
    output:
        "results/kallisto/{sample}/abundance.tsv"
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{sample}"
    threads:
        1
    conda:
        "../envs/kallisto.yaml"
    log:
        "results/logs/kallisto/{sample}.log"
    message:
        "Perform pseudoalignment and quantify transcript abundance for {wildcards.sample}."
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"

rule tx2gene:
    input:
        gtf = config["path"]["genome_gtf"]
    output:
        tsv = "data/references/tx2gene.tsv"
    conda:
        "../envs/python.yaml"
    message:
        "Convert transcript IDs to gene IDs."
    script:
        "../scripts/tx2gene.py"

rule merge_kallisto_quant:
    input:
        quant = expand(rules.kallisto_quant.output, sample=SAMPLES),
        tx2gene = rules.tx2gene.output.tsv
    output:
        tpm = "results/kallisto/tpm.tsv"
    conda:
        "../envs/python.yaml"
    log:
        "results/logs/kallisto/merge_kallisto_quant.log"
    message: 
        "Merge kallisto quantification results into one dataframe for further analysis."
    script:
        "../scripts/merge_kallisto_quant.py"

rule pca:
    input:
        tpm = rules.merge_kallisto_quant.output.tpm
    output:
        plot = "results/pca.svg"
    params:
        design = 'data/design.tsv'
    conda:
        "../envs/python_plots.yaml"
    log:
        "results/logs/pca.log"
    message:
        "Generate a PCA plot to observe variance between samples."
    script:
        "../scripts/pca.py"

rule multiqc:
    input:
        star = expand(rules.primary_alignments.output, sample=SAMPLES),
        tpm = rules.merge_kallisto_quant.output.tpm
    output:
        html = "results/multiqc/multiqc_report.html"
    params:
        scan_dir = "results/fastqc results/star/*/*Log.final.out results/kallisto",
        outdir = directory('results/multiqc')
    log:
        "results/logs/multiqc.log"
    message:
        "Summarize analysis results for multiple tools and samples in a single report."
    conda:
        "../envs/fastqc.yaml"
    shell:
        "multiqc {params.scan_dir} "
        "--outdir {params.outdir} "
        "&> {log}"