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

rule primary_alignments:
    # More robust results with multimapped reads
    input:
        rules.star_align.output.aln
    output:
        "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.primary.bam"
    log:
        "results/logs/star/{sample}_primary.log"
    conda:
        "../envs/genomecov.yaml"
    message:
        "Keep primary alignments only for {wildcards.sample}."
    shell:
        "samtools view -b -F 256 -o {output} {input} "
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
    # TmpScale counts the number of reads in the BAM and multiplies it by 1mio, rounded to 6 digits
    shell:
        "TmpScale=$(bc <<< \"scale=6;1000000/$(samtools view -f 0 -c {input})\") && "
        "bedtools genomecov -bg -split -ibam {input} -scale $TmpScale | sort -k1,1 -k2,2n > {output}"