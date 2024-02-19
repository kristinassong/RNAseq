rule star_index:
    input:
        fasta = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        "results/STAR/index/chrNameLength.txt"
    params:
        idx = "results/STAR/index"
    threads: 
        32
    conda:
        "../envs/STAR.yaml"
    log:
        "results/logs/STAR/index.log"
    message:
        "Generate genome indexes files using STAR."
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.idx} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99 "
        "&> {log}"


rule star_align:
    input:
        fq1 = rules.trimmomatic.output.r1,
        fq2 = rules.trimmomatic.output.r2,
        chrNameLength = rules.star_index.output
    output:
        bam = "results/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    log:
        "results/logs/STAR/{sample}_align.log"
    conda:
        "../envs/STAR.yaml"
    params:
        out_prefix = "results/STAR/{sample}/{sample}_",
        idx = rules.star_index.params.idx
    threads:
        32
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.idx} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--runThreadN {threads} "
        "--readFilesCommand zcat --outReadsUnmapped Fastx "
        "--outFilterType BySJout --outStd Log --outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.out_prefix} "
        "--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 "
        "--alignEndsProtrude 5 ConcordantPair "
        "&> {log}"


rule primary_alignments:
    input:
        rules.star_align.output.bam
    output:
        "results/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.primary.bam"
    log:
        "results/logs/STAR/{sample}_primary.log"
    conda:
        "../envs/STAR.yaml"
    message:
        "Keep primary alignments only for {wildcards.sample}."
    shell:
        "samtools view -b -F 256 -o {output} {input} "
        "&> {log}"


rule picard:
    input:
        rules.primary_alignments.output
    output:
        txt = "results/picard/{sample}.isize.txt",
        pdf = "results/picard/{sample}.isize.pdf"
    log:
        "results/logs/picard/{sample}.log"
    conda:
        "../envs/picard.yaml"
    message:
        "Collect insert size distribution metrics for validating library construction for {wildcards.sample}."
    shell:
        "picard CollectInsertSizeMetrics "
        "--INPUT {input} "
        "--OUTPUT {output.txt} "
        "--Histogram_FILE {output.pdf} "
        "--MINIMUM_PCT 0.5 "
        "&> {log}"


rule genomecov:
    input:
        rules.primary_alignments.output
    output:
        "results/genomecov/{sample}.bedgraph"
    conda:
        "../envs/genomecov.yaml"
    message:
        "Report {wildcards.sample} genome coverage in BEDGRAPH format."
    shell:
        "bedtools genomecov -bg -split -ibam {input} | sort -k1,1 -k2,2n > {output}"