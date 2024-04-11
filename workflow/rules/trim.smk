rule fastqc_pretrim:
    input:
        "resources/samples/{sample}_{read}.fastq.gz"
    output:
        "results/fastqc/pretrim/{sample}_{read}_fastqc.html"
    params:
        "results/fastqc/pretrim"
    log:
        "results/logs/fastqc/pretrim/{sample}_{read}.log"
    threads:
        32
    message:
        "Quality control check on raw sequence data of {wildcards.sample}_{wildcards.read}."
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "-t {threads} "
        "{input} "
        "&> {log}"


rule trimmomatic:
    input:
        r1 = "resources/samples/{sample}_R1.fastq.gz",
        r2 = "resources/samples/{sample}_R2.fastq.gz",
    output:
        r1 = "results/trimmomatic/{sample}_R1.fastq.gz",
        r2 = "results/trimmomatic/{sample}_R2.fastq.gz",
        unpaired_r1 = "results/trimmomatic/{sample}_R1.unpaired.fastq.gz",
        unpaired_r2 = "results/trimmomatic/{sample}_R2.unpaired.fastq.gz"
    log:
        "results/logs/trimmomatic/{sample}.log"
    params:
        trimmer = "ILLUMINACLIP:resources/Adapters-PE_NextSeq.fa:2:12:10:8:true TRAILING:30 LEADING:30 MINLEN:20",
    threads:
        8
    message:
        "Trim poor quality reads in {wildcards.sample} using Trimmomatic."
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.r1} {input.r2} "
        "{output.r1} {output.unpaired_r1} "
        "{output.r2} {output.unpaired_r2} "
        "{params.trimmer} "
        "2> {log}"


rule fastqc_posttrim:
    input:
        trimmed = [rules.trimmomatic.output.r1,rules.trimmomatic.output.r2],
        r = "results/trimmomatic/{sample}_{read}.fastq.gz"
    output:
        html = "results/fastqc/posttrim/{sample}_{read}_fastqc.html"
    params:
        "results/fastqc/posttrim"
    threads:
        32
    log:
        "results/logs/fastqc/posttrim/{sample}_{read}.log"
    message:
        "Quality control check on trimmed sequence data of {wildcards.sample}_{wildcards.read}."
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "-t {threads} "
        "{input.r} "
        "&> {log}"