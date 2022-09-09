conda activate cyclomics-nf-cyclomicseq
conda env export > environment.yml
docker build . -t damicyclomics/cyclomicseq:0.7.0
docker tag damicyclomics/cyclomicseq:0.7.0 damicyclomics/cyclomicseq:latest
docker push damicyclomics/cyclomicseq:0.7.0
docker push damicyclomics/cyclomicseq:latest