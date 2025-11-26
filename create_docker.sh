docker build . -t cyclomics/cyclomicsseq:2.0.0-rc1
docker tag cyclomics/cyclomicsseq:2.0.0-rc1 cyclomics/cyclomicsseq:latest
docker push cyclomics/cyclomicsseq:2.0.0-rc1
docker push cyclomics/cyclomicsseq:latest
