rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    wrapper:
        "v1.14.1/bio/star/index"

rule star_align:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1="reads/{sample}_R1.1.fastq",
        fq2="reads/{sample}_R2.1.fastq",  #optional
        idx="resources/star_genome",
    output:
        # see STAR manual for additional output files
        aln="star/{sample}/Aligned.out.bam",
        reads_per_gene="star/{sample}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}.log",
    params:
        # specific parameters to work well with arriba
        extra="--quantMode GeneCounts --sjdbGTFfile resources/genome.gtf"
        " --outSAMtype BAM Unsorted --chimSegmentMin 10 --chimOutType WithinBAM SoftClip"
        " --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0"
        " --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3",
    threads: 12
    wrapper:
        "v1.14.1/bio/star/align"