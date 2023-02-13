#!/bin/bash

### NEED TO RUN COMMANDS INDIVIDUALLY FOR NOW ###

run_dir=$SCRATCH/RNAseq/results/genomecov
chrom_sizes=$run_dir/hg38.chrom.sizes

# Get actual chrom.sizes information from UCSC
fetchChromSizes hg38 > $chrom_sizes
# change 'chrM' to 'chrMT'
sed -i 's/chrM/chrMT/' $chrom_sizes

# For all bedgraph files in genomecov directory
for f in $run_dir/*.bedgraph
do
    fbname=$(basename "$f" .bedgraph)
    f_chr=$run_dir/$fbname.chr.bedgraph
    bwname=$run_dir/$fbname.bw
    # add 'chr' to first column 
    sed 's/^/chr/' $f > $f_chr &&
done

# Convert bedgraph to bigwig format
for f in $run_dir/*.chr.bedgraph
do
    fbname=$(basename "$f" .bedgraph)
    bedGraphToBigWig $f_chr $chrom_sizes $bwname
done

# Remove intermediate files
rm $run_dir/*.chr.bedgraph