rule star_index:
    input:
        fasta = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        directory("data/references/star_index")
    threads: 
        8
    conda:
        "../envs/star.yaml"
    log:
        "results/logs/star/index.log"
    message:
        "Generate genome indexes files using STAR."
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99 "
        "&> {log}"

rule star_align:
    input:
        fq1 = rules.trimmomatic.output.r1,
        fq2 = rules.trimmomatic.output.r2,
        idx = rules.star_index.output
    output:
        aln = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    log:
        "results/logs/star/{sample}.log"
    conda:
        "../envs/star.yaml"
    params:
        extra = config["params"]["star"],
        out_prefix = "results/star/{sample}/{sample}_"
    threads:
        8
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--runThreadN {threads} "
        "--readFilesCommand zcat --outReadsUnmapped Fastx "
        "--outFilterType BySJout --outStd Log --outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.out_prefix} "
        "{params.extra} "
        "&> {log}"
"""
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
"""