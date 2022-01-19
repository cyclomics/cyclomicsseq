# CycloSeq

Pipeline to process Cyclomics data.

## Dependencies

Nextflow

## Usage

```bash
nextflow cyclomics/cycloseq -v v0.0.1dev
```

## Informed pipeline

There is an informed pipeline that takes as much prior information into account as possible. This call be called using:

```bash
nextflow run workflows/informed.nf --input_read_dir '/some/path' --read_pattern 'fastq_pass/*.fastq' --output_dir 'testing' -resume
```
