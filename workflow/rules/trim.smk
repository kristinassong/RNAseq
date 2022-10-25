rule trimmomatic:
    input:
        r1 = os.path.join(config["path"]["fastq"],"{sample}_R1_001.fastq.gz"),
        r2 = os.path.join(config["path"]["fastq"],"{sample}_R2_001.fastq.gz")
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
    wrapper:
        "v1.17.4/bio/trimmomatic/pe"
    message:
        "Filter poor quality reads in {sample} using Trimmomatic."