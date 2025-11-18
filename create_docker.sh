docker build . -t cyclomics/cyclomicsseq:0.9.2
docker tag cyclomics/cyclomicsseq:0.9.2 cyclomics/cyclomicsseq:latest
docker push cyclomics/cyclomicsseq:0.9.2
docker push cyclomics/cyclomicsseq:latest
