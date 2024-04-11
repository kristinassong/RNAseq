# RNA-seq Analysis
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) 🐍 pipeline to analyze RNA-seq expression data.

Author: [Kristina Sungeun Song](mailto:kristina.song@usherbrooke.ca)

## Installing Snakemake
Please follow the [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install the Snakemake workflow management tool. We recommend using `Conda/Mamba` to install Snakemake.

This Snakemake workflow has been tested with `v7.32.4`.

## Running the Snakemake workflow
For a dry-run of this Snakemake workflow, simply run the following code from `RNAseq/`.
```
snakemake -n
```
To run this Snakemake workflow, simply run the following code from `RNAseq/`.
```
snakemake --profile profile_slurm

OR

snakemake --profile profile_local
```

### Steps
1. Trimming: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. Read mapping: [Kallisto](https://pachterlab.github.io/kallisto/), [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)
3. Read quality check: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), Picard, PCA, [MultiQC](https://multiqc.info/docs/)
4. Differential expression analysis: [DESeq2](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
5. Alternative splicing analysis: [rMATS](https://rnaseq-mats.sourceforge.io/) (paired & non-paired)