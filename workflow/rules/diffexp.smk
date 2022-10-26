rule deseq2_init:
    input:
    output:
    log:
        "results/logs/deseq2/initialization.log"
    script:
        "workflow/scripts/deseq2-init.R"

rule deseq2:
    input:
    output:
    log:
        "results/logs/kallisto/{?????????}.log"
    script:
        "workflow/scripts/deseq2.R"