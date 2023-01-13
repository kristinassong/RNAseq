rule majiq_build:
    input:
        gff3 = config["path"]["gff3"],
        bai = expand(rules.bam_index.output, sample=SAMPLES)
    output:
        splicegraph = "results/majiq/build/splicegraph.sql",
        majiq_files = expand('results/majiq/build/{sample}_Aligned.sortedByCoord.out.primary.majiq',sample=SAMPLES)
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

rule majiq_psi_quant:
    input:
        build_dir = rules.majiq_build.output.majiq_files
    output:
        voila = "results/majiq/psi_quant/{cond}.psi.voila"
    params:
        outdir = directory("results/majiq/psi_quant"),
        group = "--name {cond}",
        majiq_files = get_majiq_files_psi,
        majiq_tool = config["tools"]["majiq_voila"]
    log:
        "results/logs/majiq/psi_quant_{cond}.log"
    message:
        "Quantify LSV candidates given by the MAJIQ builder for the {wildcards.cond} group."
    shell:
        "source {params.majiq_tool} && "
        "majiq psi {params.majiq_files} -j 8 --output {params.outdir} {params.group} "
        "&> {log} && deactivate"

rule majiq_deltapsi_quant:
    input:
        build_dir = rules.majiq_build.output.majiq_files
    output:
        voila = "results/majiq/deltapsi_quant/{comp}.deltapsi.voila"
    params:
        outdir = directory("results/majiq/deltapsi_quant"),
        group = majiq_deltapsi_name_format,
        majiq_tool = config["tools"]["majiq_voila"],
        grp1_majiq_files = majiq_cond1,
        grp2_majiq_files = majiq_cond2
    log:
        "results/logs/majiq/deltapsi_quant_{comp}.log"
    message:
        "Quantify differential splicing between two different groups: {wildcards.comp}."
    shell:
        "source {params.majiq_tool} && "
        "majiq deltapsi -grp1 {params.grp1_majiq_files} -grp2 {params.grp2_majiq_files} "
        "-j 8 -o {params.outdir} --name {params.group} "
        "&> {log} && deactivate"

rule voila_tsv_psi:
    input:
        voila = rules.majiq_psi_quant.output.voila,
        splicegraph = rules.majiq_build.output.splicegraph
    output:
        tsv = "results/voila/psi/{cond}.psi.tsv" 
    params:
        majiq_tool = config["tools"]["majiq_voila"]
    log:
        "results/logs/voila/psi_{cond}.log"
    message:
        "Generate VOILA tsv file for {wildcards.cond} PSI computation."
    shell:
        "source {params.majiq_tool} && "
        "voila tsv {input.splicegraph} {input.voila} -f {output.tsv} "
        "&> {log} && deactivate"
    
rule voila_tsv_deltapsi:
    input:
        splicegraph = rules.majiq_build.output.splicegraph,
        voila = rules.majiq_deltapsi_quant.output.voila
    output:
        tsv = "results/voila/deltapsi/{comp}.deltapsi.tsv"
    params:
        majiq_tool = config["tools"]["majiq_voila"]
    log:
        "results/logs/voila/deltapsi_{comp}.log"
    message:
        "Generate VOILA tsv file for {wildcards.comp} delta PSI computation."
    shell:
        "source {params.majiq_tool} && "
        "voila tsv {input.splicegraph} {input.voila} -f {output.tsv} "
        "&> {log} && deactivate"

rule tsv_psi_filtered:
    input:
        tsv = rules.voila_tsv_psi.output.tsv,
        exp_genes = rules.merge_kallisto_quant.output.tpm
    output:
        filtered_tsv = "results/voila/psi/{cond}.psi.filtered.tsv" 
    log:
        "results/logs/voila/psi_{cond}_filtered.log"
    message:
        "Keep VOILA psi results only for the genes that are expressed."
    script:
        "../scripts/voila_filter.py"

rule tsv_deltapsi_filtered:
    input:
        tsv = rules.voila_tsv_deltapsi.output.tsv,
        exp_genes = rules.merge_kallisto_quant.output.tpm
    output:
        filtered_tsv = "results/voila/deltapsi/{comp}.deltapsi.filtered.tsv" 
    log:
        "results/logs/voila/deltapsi_{comp}_filtered.log"
    message:
        "Keep VOILA deltapsi results only for the genes that are expressed."
    script:
        "../scripts/voila_filter.py"
