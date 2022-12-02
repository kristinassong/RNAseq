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
        bam = expand(rules.primary_alignments.output, sample=SAMPLES),
        bai = expand(rules.bam_index.output, sample=SAMPLES)
    output:
        directory("results/majiq/build")
    params:
        "data/majiq.conf"
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/majiq/build.log"
    message:
        "Analyze RNA-seq data to detect LSV candidates."
    shell:
        "majiq build {input.gff3} -c {params} -j 8 -o {output} "
        "&> {log}"

rule majiq_psi_quant:
    input:
        rules.majiq_build.output,
        get_majiq_psi_quant_inputs   
    output:
        expand("results/majiq/psi_quant/{cond}.psi.pickle",cond=condition)
    params:
        directory("results/majiq/psi_quant")
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/majiq/"
    message:
        "Quantify LSV candidates given by MAJIQ builder."
    shell:
        "majiq psi {input} --nthreads 8 --output {params} --name"
        "&> {log}"

rule voila_psi:
    input:
    output:
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/voila"
    message:
        "Generate interactive summaries to display MAJIQ computations and quantifications."
    shell:
        "voila psi {output} --genes-files -o "
        "&> {log}"

rule deltapsi_quant:
    input:
    output:
    conda:
        "../envs/majiq_voila.yaml"
    log:
        "results/logs/majiq/"
    message:
        "Quantify differential splicing between two different groups."
    shell:
        "majiq deltapsi -grp1 -grp2 -j 8 -o {output} -n "

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