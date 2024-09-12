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
        #summary = rules.rmats.output.summary,
        summary = 'results/rmats/SNORD22_OVE-pcDNA/raw/summary.txt',
        tpm = rules.merge_kallisto_quant.output.tpm
    output:
        result = 'results/rmats/{comp}/filtered/SE.tsv'
    params:
        dir = directory("results/rmats/{comp}"),
        gtf = config["path"]["genome_gtf"],
        fdr = 0.1, # recommended threshold: <=1%
        deltapsi = 0.05, # recommended threshold: >=5%
        rc = 10, # recommended thershold: 10
        basepsi_low = 0.05, # recommended thershold: 0.05
        basepsi_high = 0.95 # recommended thershold: 0.95
    conda:
        "../envs/python.yaml"
    message:
        "Filter raw rMATS output for {wildcards.comp}."
    script:
        "../scripts/filter_rmats.py"


rule rmats_paired_env:
    # REQUIREMENTS:
    # 1) git clone https://github.com/Xinglab/rmats-turbo.git to resources/
    # 2) git clone https://github.com/Xinglab/PAIRADISE.git to resources/rmats-turbo
    output:
        "results/logs/rmats_paired_env.log"
    conda:
        "../envs/rmats.yaml"
    message:
        "Install rMATS dependencies for --paired-stats."
    shell:
        "cd resources/rmats-turbo && "
        "Rscript install_r_deps.R paired && "
        "echo \'Installation complete.\' > ../../{output}"


rule rmats_paired:
    input:
        bam = expand(rules.primary_alignments.output, sample=SAMPLES),
        group1 = "resources/rmats_b1.txt",
        group2 = "resources/rmats_b2.txt",
        gtf = config["path"]["genome_gtf"],
        log = rules.rmats_paired_env.output
    output:
        outdir = directory("results/rmats_paired/{comp}/raw"),
        tmpdir = directory("results/rmats_paired/{comp}/tmp"),
        summary = "results/rmats_paired/{comp}/raw/summary.txt"
    params:
        readlength = config["params"]["readlength"]
    conda:
        "../envs/rmats.yaml"
    message:
        "Run rMATS for {wildcards.comp} using --paired-stats."
    shell:
        "rmats.py --b1 {input.group1} --b2 {input.group2} "
        "--gtf {input.gtf} -t paired --readLength {params.readlength} --variable-read-length --paired-stats "
        "--nthread 4 --od {output.outdir} --tmp {output.tmpdir}"


rule filter_rmats_paired:
    input:
        #summary = rules.rmats_paired.output.summary,
        summary = 'results/rmats_paired/SNORD22_OVE-pcDNA/raw/summary.txt',
        tpm = rules.merge_kallisto_quant.output.tpm
    output:
        result = 'results/rmats_paired/{comp}/filtered/SE.tsv'
    params:
        dir = directory("results/rmats_paired/{comp}"),
        gtf = config["path"]["genome_gtf"],
        fdr = 0.1, # recommended threshold: <=1%
        deltapsi = 0.05, # recommended threshold: >=5%
        rc = 10, # recommended thershold: 10
        basepsi_low = 0.05, # recommended thershold: 0.05
        basepsi_high = 0.95 # recommended thershold: 0.95
    conda:
        "../envs/python.yaml"
    message:
        "Filter raw rMATS --paired-stats output for {wildcards.comp}."
    script:
        "../scripts/filter_rmats.py"