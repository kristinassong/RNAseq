rule rmats:
    input:
        bam = expand(rules.primary_alignments.output, sample=SAMPLES),
        group1 = "resources/rmats_b1.txt",
        group2 = "resources/rmats_b2.txt",
        gtf = config["path"]["genome_gtf"]
    output:
        outdir = directory("results/rmats/{comp}/raw"),
        tmpdir = directory("results/rmats/{comp}/tmp"),
        summary = "results/rmats/{comp}/raw/summary.txt"
    params:
        readlength = config["params"]["readlength"]
    conda:
        "../envs/rmats.yaml"
    message:
        "Run rMATS for {wildcards.comp}."
    shell:
        "rmats.py --b1 {input.group1} --b2 {input.group2} "
        "--gtf {input.gtf} -t paired --readLength {params.readlength} --variable-read-length "
        "--nthread 4 --od {output.outdir} --tmp {output.tmpdir}"


rule filter_rmats:
    input:
        summary = rules.rmats.output.summary,
        tpm = rules.merge_kallisto_quant.output.tpm
    output:
        result = 'results/rmats/{comp}/filtered/SE.tsv'
    params:
        dir = directory("results/rmats/{comp}"),
        gtf = config["path"]["genome_gtf"],
        fdr = 0.05,
        deltapsi = 0.10
    conda:
        "../envs/python.yaml"
    message:
        "Filter raw rMATS output for {wildcards.comp}."
    script:
        "../scripts/filter_rmats.py"
