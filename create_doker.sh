conda activate cyclomics-nf-cyclomicseq
conda env export > environment.yml
docker build . -t damicyclomics/cyclomicseq:0.7.1
docker tag damicyclomics/cyclomicseq:0.7.1 damicyclomics/cyclomicseq:latest
docker push damicyclomics/cyclomicseq:0.7.1
docker push damicyclomics/cyclomicseq:latest