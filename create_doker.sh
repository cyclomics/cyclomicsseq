conda activate cyclomics-nf-cyclomicseq
conda env export > environment.yml
docker build . -t damicyclomics/cyclomicseq:0.0.1
