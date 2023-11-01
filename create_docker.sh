docker build . -t damicyclomics/cyclomicseq:0.8.0
docker tag damicyclomics/cyclomicseq:0.8.0 damicyclomics/cyclomicseq:latest
docker push damicyclomics/cyclomicseq:0.8.0
docker push damicyclomics/cyclomicseq:latest
