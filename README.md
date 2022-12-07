# RNAseq
A snakemake pipeline to analyze RNA-seq expression data

### Preinstallation of HTSlib library for MAJIQ and VOILA
```
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
./configure --prefix=/home/kris98/scratch/RNAseq/workflow/envs/htslib-1.9/htslib
make
make install
$ export HTSLIB_LIBRARY_DIR=/home/kris98/scratch/RNAseq/workflow/envs/htslib-1.9/htslib/lib
$ export HTSLIB_INCLUDE_DIR=/home/kris98/scratch/RNAseq/workflow/envs/htslib-1.9/htslib/include
```

### Install MAJIQ using venv
```
python3.8 -m venv majiq_voila
source majiq_voila/bin/activate
pip install git+https://bitbucket.org/biociphers/majiq_academic.git"
deactivate
```
