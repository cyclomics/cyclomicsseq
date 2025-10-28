docker build . -t cyclomics/cyclomicsseq:0.9.1
docker tag cyclomics/cyclomicsseq:0.9.1 cyclomics/cyclomicsseq:latest
docker push cyclomics/cyclomicsseq:0.9.1
docker push cyclomics/cyclomicsseq:latest
