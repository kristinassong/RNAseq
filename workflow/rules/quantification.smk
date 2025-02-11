rule build_transcriptome:
    input:
        genome = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        'resources/transcriptome.fa'
    log:
        "results/logs/gffread.log"
    conda:
        "../envs/gffread.yaml"
    message:
        "Build a reference transcriptome using gffread."
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output} &> {log}"


rule kallisto_index:
    input:
        rules.build_transcriptome.output
    output:
        "results/kallisto/kallisto.idx"
    params:
        31
    conda:
        "../envs/kallisto.yaml"
    log:
        "results/logs/kallisto/index.log"
    message:
        "Build a Kallisto index from the transcriptome FASTA file."
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
        tsv = "resources/tx2gene.tsv"
    conda:
        "../envs/python.yaml"
    message:
        "Convert transcript IDs to gene IDs."
    script:
        "../scripts/tx2gene.py"


rule merge_kallisto_quant:
    input:
        quant = expand(rules.kallisto_quant.output, sample=SAMPLES),
        tx2gene = rules.tx2gene.output.tsv,
        gtf = config["path"]["genome_gtf"]
    output:
        tpm = "results/kallisto/tpm_{comp}.tsv",
        counts = "results/kallisto/counts_{comp}.tsv"
    conda:
        "../envs/python.yaml"
    message: 
        "Merge kallisto quantification results into one dataframe for further analysis."
    script:
        "../scripts/merge_kallisto_quant.py"