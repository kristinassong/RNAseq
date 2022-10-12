rule kallisto_index:
    input:
        fasta="{transcriptome}.fasta",
    output:
        index="{transcriptome}.idx",
    log:
        "results/logs/kallisto/index_{transcriptome}.log",
    threads: 1
    message:
        "Builds an index from {input.fasta}."
    wrapper:
        "v1.14.1/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq=["reads/{exp}_R1.fastq", "reads/{exp}_R2.fastq"],
        index="index/transcriptome.idx",
    output:
        directory("quant_results_{exp}"),
    log:
        "results/logs/kallisto/quant_{exp}.log",
    threads: 1
    wrapper:
        "v1.14.1/bio/kallisto/quant"