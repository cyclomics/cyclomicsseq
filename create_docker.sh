docker build . -t cyclomics/cyclomicsseq:0.9.0
docker tag cyclomics/cyclomicsseq:0.9.0 cyclomics/cyclomicsseq:latest
docker push cyclomics/cyclomicsseq:0.9.0
docker push cyclomics/cyclomicsseq:latest
