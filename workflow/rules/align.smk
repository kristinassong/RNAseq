rule star_index:
    input:
        fasta = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        directory(config["path"]["star_index"])
    threads: 
        8
    params:
        extra = "--sjdbOverhang 99"
    log:
        "results/logs/star/index.log"
    message:
        "Generate genome indexes files using STAR."
    wrapper:
        "v1.17.4/bio/star/index"

rule star_align:
    input:
        fq1 = rules.trimmomatic.output.r1,
        fq2 = rules.trimmomatic.output.r2,
        idx = rules.star_index.output
    output:
        aln = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log = "results/star/{sample}/Log.out"
    log:
        "results/logs/star/{sample}.log"
    params:
        extra = config["params"]["star"]
    threads:
        8
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    wrapper:
        "v1.17.4/bio/star/align"

rule genomecov:
    input:
        rules.star_align.output.aln
    output:
        temp("results/genomecov/{sample}.temp.bedgraph")
    log:
        "results/logs/genomecov/{sample}.log"
    params:
        "-bg -split"
    message:
        "Report {wildcards.sample} genome coverage in BEDGRAPH format."
    wrapper:
        "v1.17.4/bio/bedtools/genomecov"

rule genomecov_sorted:
    input:
        rules.genomecov.output
    output:
        "results/genomecov/{sample}.bedgraph"
    message:
        "Sort {wildcards.sample} genome coverage bedgraph by chromosome and start position."
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"