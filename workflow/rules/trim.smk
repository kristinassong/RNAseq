rule trimmomatic:
    input:
    output:
    logs:
        "results/logs/trimmomatic/{FIX}.log"
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        