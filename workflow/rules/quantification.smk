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