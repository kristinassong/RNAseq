rule bcl2fastq2:
    input:
        config['samples']
    output:
        fastq = directory(Path("results/bcl2fastq2", '{id}'))
        "results/bcl2fastq2/"
    conda:
        "envs/bcl2fastq2.yaml"
    log:
        "results/logs/bcl2fastq2/{id}.log"
    message:
        "Convert bcl files to fastq files."
    shell:
        "bcl2fastq2 -R {FIXX} -o {output} --sample-sheet {input} "
        "--minimum-trimmed-read-length 13 --mask-short-adapter-reads 13 --no-lane-splitting -p 4"