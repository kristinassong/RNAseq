rule fastqc_pretrim:
    input:
        os.path.join(config["path"]["fastq"],"{sample}_{read}_001.fastq.gz")
    output:
        html = "results/fastqc/pretrim/{sample}_{read}.html",
        zip = "results/fastqc/pretrim/{sample}_{read}_fastqc.zip"
    log:
        "results/logs/fastqc/pretrim/{sample}_{read}.log"
    message:
        "Quality control check on raw sequence data of {wildcards.sample}_{wildcards.read}."
    wrapper:
        "v1.17.4/bio/fastqc"

rule fastqc_posttrim_R1:
    input:
        rules.trimmomatic.output.r1
    output:
        html = "results/fastqc/posttrim/{sample}_R1.html",
        zip = "results/fastqc/posttrim/{sample}_R1_fastqc.zip"
    log:
        "results/logs/fastqc/posttrim/{sample}_R1.log"
    message:
        "Quality control check on trimmed sequence data of {wildcards.sample}_R1."
    wrapper:
        "v1.17.4/bio/fastqc"

rule fastqc_posttrim_R2:
    input:
        rules.trimmomatic.output.r2
    output:
        html = "results/fastqc/posttrim/{sample}_R2.html",
        zip = "results/fastqc/posttrim/{sample}_R2_fastqc.zip"
    log:
        "results/logs/fastqc/posttrim/{sample}_R2.log"
    message:
        "Quality control check on trimmed sequence data of {wildcards.sample}_R2."
    wrapper:
        "v1.17.4/bio/fastqc"