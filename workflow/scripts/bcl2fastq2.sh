#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.song@usherbrooke.ca
#SBATCH --time=5:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --error=%u.%x-%A[%a].err
#SBATCH --mem=48000M
#SBATCH --cpus-per-task=12

### Modified from Gabrielle's script 

run=../../data/samples/220923_NB502083_0194_AHVGMKBGXM
path=$run/bcl
sample_sheet=$path/sample_sheet.csv
outpath=$run/fastq

mkdir -p $outpath &&

bcl2fastq -R $path -o $outpath --sample-sheet $sample_sheet \
--minimum-trimmed-read-length 13 --mask-short-adapter-reads 13 \
--no-lane-splitting -p 4 
