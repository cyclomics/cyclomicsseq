docker build . -t damicyclomics/cyclomicseq:0.7.2
docker tag damicyclomics/cyclomicseq:0.7.2 damicyclomics/cyclomicseq:latest
docker push damicyclomics/cyclomicseq:0.7.2
docker push damicyclomics/cyclomicseq:latest
