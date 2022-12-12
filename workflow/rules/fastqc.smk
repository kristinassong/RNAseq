rule fastqc_pretrim:
    input:
        os.path.join(config["path"]["fastq"],"{sample}_{read}_001.fastq.gz")
    output:
        html = "results/fastqc/pretrim/{sample}_{read}_001_fastqc.html",
        zip = "results/fastqc/pretrim/{sample}_{read}_001_fastqc.zip"
    params:
        "results/fastqc/pretrim"
    log:
        "results/logs/fastqc/pretrim/{sample}_{read}.log"
    message:
        "Quality control check on raw sequence data of {wildcards.sample}_{wildcards.read}."
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input} "
        "&> {log}"

rule pretrim_fname:
    input:
        html = rules.fastqc_pretrim.output.html,
        zip = rules.fastqc_pretrim.output.zip
    output:
        html = "results/fastqc/pretrim/{sample}_{read}_fastqc.html",
        zip = "results/fastqc/pretrim/{sample}_{read}_fastqc.zip"
    shell:
        "touch {output.html} && touch {output.zip} && "
        "mv {input.html} {output.html} && "
        "mv {input.zip} {output.zip}"

rule fastqc_posttrim_R1:
    input:
        rules.trimmomatic.output.r1
    output:
        html = "results/fastqc/posttrim/{sample}_R1_fastqc.html",
        zip = "results/fastqc/posttrim/{sample}_R1_fastqc.zip"
    params:
        "results/fastqc/posttrim"
    log:
        "results/logs/fastqc/posttrim/{sample}_R1.log"
    message:
        "Quality control check on trimmed sequence data of {wildcards.sample}_R1."
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input} "
        "&> {log}"

rule fastqc_posttrim_R2:
    input:
        rules.trimmomatic.output.r2
    output:
        html = "results/fastqc/posttrim/{sample}_R2_fastqc.html",
        zip = "results/fastqc/posttrim/{sample}_R2_fastqc.zip"
    params:
        "results/fastqc/posttrim"
    log:
        "results/logs/fastqc/posttrim/{sample}_R2.log"
    message:
        "Quality control check on trimmed sequence data of {wildcards.sample}_R2."
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input} "
        "&> {log}"
