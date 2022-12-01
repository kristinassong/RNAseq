rule picard:
    input:
        rules.primary_alignments.output
    output:
        txt = "results/picard/{sample}.isize.txt",
        pdf = "results/picard/{sample}.isize.pdf"
    log:
        "results/logs/picard/{sample}.log"
    params:
        config["params"]["picard"]
    conda:
        "../envs/picard.yaml"
    message:
        "Collect insert size distribution metrics for validating library construction for {wildcards.sample}."
    shell:
        "picard CollectInsertSizeMetrics "
        "--INPUT {input} "
        "--OUTPUT {output.txt} "
        "--Histogram_FILE {output.pdf} "
        "{params} "
        "&> {log}"