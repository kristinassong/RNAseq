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

rule voila_deltapsi_categorize_events:
    input:
        splicegraph = rules.majiq_build.output.splicegraph,
        voila = rules.majiq_deltapsi_quant.output.voila
    output:
        summary = "results/voila/deltapsi/event_types/{comp}/raw/summary.tsv"
    params:
        outdir = directory("results/voila/deltapsi/event_types/{comp}/raw"),
        majiq_tool = config["tools"]["majiq_voila"]
    log:
        "results/logs/voila/deltapsi_event_types_{comp}_raw.log"
    message:
        "Categorize splicing events within modules using VOILA for {wildcards.comp}."
    shell:
        "source {params.majiq_tool} && "
        "voila modulize {input.splicegraph} {input.voila} "
        "-d {params.outdir} -j 8 --overwrite "
        "&> {log} && deactivate"

rule filter_genes:
    input:
        summary = rules.voila_deltapsi_categorize_events.output.summary,
        exp_genes = rules.merge_kallisto_quant.output.tpm
    output:
        bar_chart = "results/voila/deltapsi/event_types/{comp}/filtered/bar_chart.svg"
    params:
        indir = directory("results/voila/deltapsi/event_types/{comp}/raw"),
        outdir = directory("results/voila/deltapsi/event_types/{comp}/filtered")
    conda:
        "../envs/python_plots.yaml"
    log:
        "results/logs/voila/deltapsi_{comp}_event_types_filtered.log"
    message:
        "Keep VOILA events only for the genes that are expressed and generate a bar chart to represent the number of events for each splicing event type."
    script:
        "../scripts/voila_filter.py"