# CycloSeq

Pipeline to process Cyclomics data.

This pipeline uses prior information from the backbone to increase the effectiveness of consensus calling from the circular DNA protocol by Cyclomics.  

## Dependencies

- Nextflow
- Docker
- Access to the damicyclomics dockerhub repo

### data requirements

- data output by ONT Guppy-hac

## Usage

### example

```bash
nextflow run cyclomics/cycloseq --input_read_dir '/some/path' --read_pattern 'fastq_pass/*.fastq' --output_dir 'testing' -resume
```

### Flag descriptions

|flag                           | info  |
|-------------------------------|---|
|--input_read_dir               | Directory where the output of Guppy is located, e.g.: "/data/guppy/exp001/".|
|--read_pattern                 | Regex pattern to look for fastq's in the read directory, defaults to: "fastq_pass/**.{fq,fastq}".|
|--sequencing_quality_summary   | Regex pattern for the summary file, default: "sequencing_summary*.txt".|
|--backbone_fasta               | Path to the fasta file containing backbone sequences.|
|--backbone_name                | Name of the sequence to extract from the backbone fasta excluding the starting ">", e.g.:"BB22". |
