rule build_transcriptome:
    input:
        genome = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        config["path"]["transcriptome"]
    conda:
        "workflow/envs/gffread.yaml"
    message:
        "Build a transcriptome using gffread."
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output}"

rule kallisto_index:
    input:
        fasta = rules.build_transcriptome.output
    output:
        index = config["path"]["kallisto_index"]
    log:
        "results/logs/kallisto/index.log"
    threads: 1
    params: 
        extra = "--kmer-size=31"
    message:
        "Builds an index from {input.fasta}."
    wrapper:
        "v1.17.4/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq = [rules.trimmomatic.output.r1, rules.trimmomatic.output.r2],
        index = rules.kallisto_index.output.index,
    output:
        directory("results/kallisto/{sample}"),
    log:
        "results/logs/kallisto/{sample}.log",
    threads: 1
    params:
        extra = "--bias --bootstrap-samples=50"
    message:
        "Quantify abundance of {wildcards.sample} reads."
    wrapper:
        "v1.17.4/bio/kallisto/quant"