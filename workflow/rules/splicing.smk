rule gtf2gff3:
    input:
        config["path"]["genome_gtf"]
    output:
        "data/references/hg38_Ensembl_V101_Scottlab_2020.gff3"
    conda:
        "../envs/gffread.yaml"
    message:
        "Convert annotation format from gtf to gff3."
    shell:
        "gffread -E {input} -o- > {output}"

rule majiq_build:
    input:
        gff3 = rules.gtf2gff3.output,
        bai = expand(rules.bam_index.output, sample=SAMPLES)
    output:
        expand("results/majiq/build/{sample}_Aligned.sortedByCoord.out.primary.majiq", sample=SAMPLES)
    params:
        config = "data/majiq.conf",
        outdir = directory("results/majiq/build"),
        majiq = config["tools"]["majiq_voila"]
    log:
        "results/logs/majiq/build.log"
    message:
        "Analyze RNA-seq data to detect LSV candidates using MAJIQ. "
        "Preinstallation of MAJIQ is required prior to running this rule."
    shell:
        "source {params.majiq} && "
        "majiq build {input.gff3} -c {params.config} -j 8 -o {params.outdir} "
        "&> {log} && deactivate"
"""
rule majiq_psi_quant:
    input:
        build_dir = rules.majiq_build.output,
        majiq_files = get_majiq_files 
    output:
        pickle = "results/majiq/psi_quant/{cond}.psi.pickle"
    params:
        outdir = directory("results/majiq/psi_quant"),
        group = "--name {cond}"
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/majiq/psi_quant_{cond}.log"
    message:
        "Quantify LSV candidates given by MAJIQ builder for {wildcards.cond} group."
    shell:
        "majiq psi {input.majiq_files} --nthreads 8 --output {params.outdir} {params.group} "
        "&> {log}"

rule voila_psi:
    input:
        pickle = rules.majiq_psi_quant.output.pickle,
        splicegraph = get_splicegraph_files
    output:
        index = "results/voila/psi/{cond}/index.html"
    params:
        directory("results/voila/psi/{cond}")
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/voila/psi_{cond}.log"
    message:
        "Generate interactive summaries to display MAJIQ computations and quantifications for {wildcards.cond}."
    shell:
        "voila psi {input.pickle} --genes-files {input.splicegraph} -o {params} "
        "&> {log}"

rule deltapsi_quant:
    input:
        build_dir = rules.majiq_build.output,
        grp1_majiq_files = get_majiq_files,
        grp2_majiq_files = get_majiq_files
    output:
        pickle = "results/majiq/deltapsi_quant/{cond1}_{cond2}.deltapsi.pickle"
    params:
        outdir = directory("results/majiq/deltapsi_quant"),
        group = "--name {cond1} {cond2}"
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/majiq/deltapsi_quant_{cond1}_{cond2}.log"
    message:
        "Quantify differential splicing between two different groups: {wildcards.cond1} and {wildcards.cond2}."
    shell:
        "majiq deltapsi -grp1 {input.grp1_majiq_files} -grp2 {input.grp2_majiq_files} "
        "-j 8 -o {params.outdir} {params.group} "
        "&> {log}"

rule voila_deltapsi:
    input:
    output:
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/voila"
    message:
        "Generate interactive summaries to display MAJIQ computations and quantifications."
    shell:
        "voila deltapsi {output} --genes-files -o "
        "&> {log}"
"""