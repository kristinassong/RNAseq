rule trimmomatic:
    input:
        r1 = os.path.join(config["path"]["fastq"],"{sample}_R1_001.fastq.gz"),
        r2 = os.path.join(config["path"]["fastq"],"{sample}_R2_001.fastq.gz"),
    output:
        r1 = "results/trimmomatic/{sample}_R1.fastq.gz",
        r2 = "results/trimmomatic/{sample}_R2.fastq.gz",
        unpaired_r1 = "results/trimmomatic/{sample}_R1.unpaired.fastq.gz",
        unpaired_r2 = "results/trimmomatic/{sample}_R2.unpaired.fastq.gz"
    log:
        "results/logs/trimmomatic/{sample}.log"
    params:
        trimmer = config["params"]["trimmomatic"],
        extra = "-phred33"
    threads:
        8
    message:
        "Filter poor quality reads in {wildcards.sample} using Trimmomatic."
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "{params.extra} "
        "{input.r1} {input.r2} "
        "{output.r1} {output.unpaired_r1} "
        "{output.r2} {output.unpaired_r2} "
        "{params.trimmer} "
        "&> {log}"