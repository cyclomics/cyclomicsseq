docker build . -t damicyclomics/cyclomicseq:0.8.1
docker tag damicyclomics/cyclomicseq:0.8.1 damicyclomics/cyclomicseq:latest
docker push damicyclomics/cyclomicseq:0.8.1
docker push damicyclomics/cyclomicseq:latest
