# CyclomicsSeq

Pipelines to process Cyclomics data.

This pipeline uses concatemeric CySeq reads, with or without backbone, as input to generate consensus reads, which will then be used to call variants over a reference.


  - [Dependencies and requirements](#dependencies)
  - [Usage](#usage)
    - [Through EPI2ME](#through-epi2me)
    - [Through command line](#through-command-line)
    - [User options](#user-options)
    - [Advanced user options](#advanced-user-options)
    - [Running on A SLURM cluster](#running-on-a-slurm-cluster)
  - [Changelog](#changelog)
  - [Dependency installation](#dependency-installation)
  - [Developer notes](#developer-notes)

## Dependencies

Click for installation instructions:

- [Nextflow](#dependency-installation) (v20.10.0 or higher)
- [Docker](#dependency-installation) or [Conda](#dependency-installation) or [Apptainer/Singularity](#dependency-installation)
- Acces to the Github repo and a valid [PAT token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)


### Data requirements

- FASTQ data output by ONT Guppy
- FASTA reference genome, ideally pre indexed by BWA to reduce runtime.


### System requirements

The pipeline expects at least 16 threads to be available and 16GB of RAM. We recommend 64 GB of RAM to decrease the runtime significantly.


### Reference genome

The pipeline has been developed with amplicons that map against the provided reference in mind.

We suggest to use the official GRCh38.p14 major release for alignment, since this works well with our variant annotation module. This release is available via the code snippet below.

``` bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

To reduce runtime pre index the reference genome with BWA, or obtain a preindexed copy from the same FTP site.


## Usage

### Through EPI2ME

This pipeline is compatible with the EPI2ME platform by ONT. Please see [ONT's installation guide](https://epi2me.nanoporetech.com/epi2me-docs/quickstart/).

Installation inside EPI2ME:
1. Go to workflows by clicking on "installed workflows", or click the workflows icon in the top bar.
2. click "Import workflow".
3. Paste "https://github.com/cyclomics/cyclomicsseq" into the text bar and click Import workflow.

Updating workflow on EPI2ME:
1.

### Through command line

In this section we assume that you have docker and nextflow installed on your system, if so running the pipeline is straightforward. You can run the pipeline directly from this repo, or pull it yourself and point nextflow towards it.

```bash
nextflow run cyclomics/cycmomicsseq -r <pipeline version> -profile docker --input_read_dir '/sequencing/20220209_1609_X3_FAS06478_0ed4361c/fastq_pass/' --output_dir '/data/myresults' --reference '/data/reference/chm13v2.fasta' --backbone BB12
```


#### Singularity

If docker is not an option, singularity (or Apptainer, as it is called since Q2 2022) is a good alternative that does not require root access and therefor used in shared compute environments.

The command becomes:

```bash
nextflow run cyclomics/cycloseq -profile singularity ...'
```

Please note that this assumes you've ran the pipeline before, if not add the -user flag as described in Usage[#Usage].

#### Conda

The pipeline is fully compatible with Conda. 
This means the full command becomes:

```bash
nextflow run cyclomics/cycloseq -profile conda ...'
```

By default it uses the environment file that is shipped with the pipeline. 
this file is located in the repo, the pipeline needs to know where this file is to run with the correct versions of the required software.


### User options

| Flag                          | Description  | Default  |
|-------------------------------|--------------|----------|
| --input_read_dir               | Directory where the output fastqs of Guppy are located, e.g.: "/data/guppy/exp001/fastq_pass". | |
| --read_pattern                 | Regex pattern to look for fastq's in the read directory. | "**.{fq,fastq,fq.gz,fastq.gz}" |
| --reference                    | Path to the reference genome to use, will ingest all index files in the same directory.| |
| --region_file                  | Path to the BED file with regions over which to run the analysis of consensus reads. | auto |
| --sample_id                    | An ID for the sample analysed. If the run is barcoded and the input directory contains many barcode subdirectories, then the barcode names will automatically be taken as sample IDs. | |
| --output_dir                   | Directory path where the results, including intermediate files, are stored. | "$HOME/Data/CyclomicsSeq" |
| --report                       | Type of report to generate at the end of the pipeline [standard, detailed, skip]. | standard |



### Advanced user options

| Flag                          | Description  | Default  |
|-------------------------------|--------------|----------|
| --max_fastq_size                | Maximum number of reads per fastq files in the analysis, only used when `--split_fastq_by_size true`. | 40000 |
| --min_repeat_count                 | Minimum number of identified repeats for reads to be considered for analysis. | 3 |
| --filtering.minimum_raw_length                  | Minimum length for reads to be considered for analysis. | 500 |
| --sequencing_summary_path      | Sequencing summary file in txt format created by MinKNOW for this experiment. | "sequencing_summary*.txt". |
| --backbone                     | Backbone used, if any [BBCS, BBCR, BB22, BB25, BB41, BB42]. | BBCS |
| --backbone_file                | File to use as backbone when --backbone is non of the available presets. eg a fasta file with a sequence with the name ">BB_custom" the name must start with BB for extraction reasons. | "$projectDir/backbones/BBCS.fasta" |
| --sequence_summary_tagging | Set whether to tag sequence summary information into BAM file annotations | false |
| --split_fastq_by_size                | Split larger fastq files into smaller ones for performance reasons. | true |
| --split_on_adapter                 | Identifies adapter sequences and splits reads when encoutered. | false |
| --include_fastq_fail                 | Include fastq_fail folder in the analysis. | false |
| --filter_by_alignment_rate        | Set whether alignments should be filtered by the rate to which they align to regions in a provided `--region_file`. | false |
| --min_align_rate                  | Minimum alignment rate of a read to designated regions of interest. Only valid if `--filter_by_alignment_rate true`. | 0.8 |
| --roi_detection.min_depth                  | Set the minimum depth required to detect a region of interest. Only valid if `--region_file auto`. | 5000 |
| --roi_detection.max_distance                 | Maximum allowed distance between detected regions for them to be merged. Only valid if `--region_file auto`. | 25 |
| --perbase.max_depth               | Set the max depth while building perbase tables. If the depth at a position is within 1% of this value, the NEAR_MAX_DEPTH output field will be set to true. | 4000000 |
| --minimap2.min_chain_score                  | Discard chains with a chaining score bewlow this value. Corresponds to minimap2 argument `-m`. | 1 |
| --minimap2.min_chain_count                  | Discard chains with number of minimizers below this value. Corresponds to minimap2 argument `-n`. | 10 |
| --minimap2.min_peak_aln_score                  |  	Minimal peak DP alignment score. Corresponds to minimap2 argument `-s`. | 20 |
| --bwamem.max_mem_occurance                 | Discard a MEM if it has more occurances than this value. Corresponds to bwa mem option `-c`. | 100 |
| --bwamem.softclip_penalty                 | Clipping penalty. Corresponds to bwa mem option `-L`. | 0.05 |
| --metadata.subsample_size                  | Maximum number of entries to take into account when plotting metadata in the report. | 10000 |
| --snp.min_ao, --indel.min_ao             | Minimum number of variant-supporting reads. | 10 |
| --snp.min_dpq, --indel.min_dpq             | Minimum positional depth after Q filtering. | 5000 |
| --snp.min_dpq_n, --indel.min_dpq_n         | Number of flanking nucleotides to the each position that will determine the window size for local maxima calculation. | 25 |
| --snp.min_dpq_ratio, --indel.min_dpq_ratio | Ratio of local depth maxima that will determine the minimum depth at each position. | 0.3 |
| --snp.min_rel_ratio, --indel.min_rel_ratio | Minimum relative ratio between forward and reverse variant-supporting reads. | 0.3 (SNP); 0.4 (Indel) |
| --snp.min_abq, --indel.min_abq             | Minimum average base quality. | 70 |
| --dynamic_vaf_params_file | Path to a YML file with dynamic parameters for variant filtering by variant allele frequency in function of total depth. | "${projectDir}/bin/variant_calling/dynamic_vaf_params.yml" |
| --max_cpus                  | Maximum number of CPUs allocated to the analysis run. | 8 |
| --max_mem_gb                  | Maximum memory in GB allocated to the analysis run. | 31 |
| --economy_mode | Set pipeline to run with a limited CPUs and memory allocation. | null |
| -profile                       | Profile with resource allocation presets and chosen environment [standard, conda, singularity, promethion]. | "standard" |


## Running on A SLURM cluster

Login to the hpc using SSH. There, start a sjob with:

```bash
srun --job-name "InteractiveJob" --cpus-per-task 16 --mem=32G --gres=tmpspace:450G --time 24:00:00 --pty bash
```

Go to the right project directory and start the pipeline as normal [through the command line](#through-command-line).


## Changelog

Please see CHANGELOG.md


## Dependency installation

### Nextflow

Download the latest version by running the example below:

```bash
wget -qO- https://get.nextflow.io | bash
```

Or see [The official Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html#installation).

### Conda

Download the latest conda version from [The official conda documentation](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

Run the below command and follow process:
```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

### Apptainer/Singularity

Download the latest version by running the example below:

```bash
wget https://github.com/apptainer/apptainer/releases/download/v1.1.0-rc.2/apptainer_1.1.0-rc.2_amd64.deb
sudo apt-get install -y ./apptainer_1.1.0-rc.2_amd64.deb
```

For the latest up to date information see their [official documentation](https://apptainer.org/docs/admin/main/installation.html)


## Developer notes

### Cycas addition to the repo

Cycas was added as a subtree using code from: https://gist.github.com/SKempin/b7857a6ff6bddb05717cc17a44091202.
This was done instead of submodule to make pulling of the repo easier for endusers and to stay compatible with `nextflow run <remote>` functionallity.


More specifically:
``` bash
git subtree add --prefix Cycas https://github.com/cyclomics/Cycas 0.4.3 --squash
```

To update, run:
``` bash
git subtree pull --prefix Cycas https://github.com/cyclomics/Cycas <tag> --squash
```
