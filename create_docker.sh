docker build . -t cyclomics/cyclomicsseq:1.2.0
docker tag cyclomics/cyclomicsseq:1.2.0 cyclomics/cyclomicsseq:latest
docker push cyclomics/cyclomicsseq:1.2.0
docker push cyclomics/cyclomicsseq:latest
