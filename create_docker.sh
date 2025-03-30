docker build . -t cyclomics/cyclomicsseq:0.8.2
docker tag cyclomics/cyclomicsseq:0.8.2 cyclomics/cyclomicsseq:latest
docker push cyclomics/cyclomicsseq:0.8.2
docker push cyclomics/cyclomicsseq:latest
