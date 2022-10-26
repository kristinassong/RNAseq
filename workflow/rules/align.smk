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
    wrapper:
        "v1.17.4/bio/star/index"
    message:
        "Generate genome indexes files using STAR."

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
    wrapper:
        "v1.17.4/bio/star/align"
    message:
        "Align {sample} reads to the reference genome using STAR."