rule fastqc_pretrim:
    input:
        os.path.join(config["path"]["fastq"],"{sample}_{R}_001.fastq.gz")
    output:
        html = "results/fastqc_pretrim/{sample}_{R}.html",
        zip = "results/fastqc_pretrim/{sample}_{R}_fastqc.zip"
    log:
        "results/logs/fastqc_pretrim/{sample}_{R}.log"
    wrapper:
        "v1.17.4/bio/fastqc"
    message:
        "Quality control check on raw sequence data of {sample}_{R}."

rule fastqc_posttrim_R1:
    input:
        rules.trimmomatic.output.r1
    output:
        html = "results/fastqc_posttrim/{sample}_R1.html",
        zip = "results/fastqc_posttrim/{sample}_R1_fastqc.zip"
    log:
        "results/logs/fastqc_posttrim/{sample}_R1.log"
    wrapper:
        "v1.17.4/bio/fastqc"
    message:
        "Quality control check on trimmed sequence data of {sample}_R1."

rule fastqc_posttrim_R2:
    input:
        rules.trimmomatic.output.r2
    output:
        html = "results/fastqc_posttrim/{sample}_R2.html",
        zip = "results/fastqc_posttrim/{sample}_R2_fastqc.zip"
    log:
        "results/logs/fastqc_posttrim/{sample}_R2.log"
    wrapper:
        "v1.17.4/bio/fastqc"
    message:
        "Quality control check on trimmed sequence data of {sample}_R2."