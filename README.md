# CycloSeq

Pipelines to process Cyclomics data.

This pipeline uses prior information from the backbone to increase the effectiveness of consensus calling from the circular DNA protocol by Cyclomics.  

## Dependencies

Click for installation instructions:

- [Nextflow](#install-dependencies)
- [Docker](#install-dependencies) or [Conda](#install-dependencies) or [Apptainer/singularity](#install-dependencies)
- Acces to the Github repo and a valid [PAT token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)


### Data requirements

- data output by ONT Guppy (SUP preferred for optimal results)
- Reference genome

## System requirements

The pipeline expects at least 16 threads to be available and 16GB of RAM. We recommend 64 GB of RAM to decrease the runtime significantly.


## Reference genome

The pipeline has been developed with amplicons that map against the provided reference in mind.

We suggest to use Grch38.p14, since this works well with the VEP that is integrated in the pipeline its available via the code snippet below.

``` bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
gunzip GCA_000001405.29_GRCh38.p14_genomic.fna.gz
```



## Usage

In this section we assume that you have docker and nextflow installed on your system, if so running the pipeline is straightforward. You can run the pipeline directly from this repo, or pull it yourself and point nextflow towards it.

If you want to run the pipeline directly from github you need to use a Personal Access Token (PAT) as the password. Click the link to see how to create a [Personal Access Token (PAT)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). Your PAT should have at least read permissions.
an example whould be:

```bash
nextflow run cyclomics/cycloseq -user <github_username> -r <current_version> -profile docker --input_read_dir '/sequencing/20220209_1609_X3_FAS06478_0ed4361c' --output_dir '/data/myresults' --reference '/data/reference/chm13v2.fasta'
```

Fill in your PAT when promted for it.


### Singularity

If docker is not an option, singularity (or Apptainer, as it is called since Q2 2022) is a good alternative that does not require root access and therefor used in shared compute environments.

The command becomes:

```bash
nextflow run cyclomics/cycloseq -profile singularity --input_read_dir '/sequencing/20220209_1609_X3_FAS06478_0ed4361c' --output_dir '/data/myresults' --reference '/data/reference/chm13v2.fasta'
```

Please note that this assumes you've ran the pipeline before, if not add the -user flag as described in Usage[#Usage].

### Conda

The pipeline is fully compatible with Conda. 
This means the full command becomes:

```bash
nextflow run cyclomics/cycloseq -profile conda --input_read_dir '/sequencing/20220209_1609_X3_FAS06478_0ed4361c' --output_dir '/data/myresults' --reference '/data/reference/chm13v2.fasta'
```

By default it uses the environment file that is shipped with the pipeline. 
this file is located in the repo, the pipeline needs to know where this file is to run with the correct versions of the required software.


### Flag descriptions

| flag                          | info  |
|-------------------------------|-------|
|--input_read_dir               | Directory where the output of Guppy is located, e.g.: "/data/guppy/exp001/".|
|--read_pattern                 | Regex pattern to look for fastq's in the read directory, defaults to: "{pass,fastq_pass}/**.{fq,fastq,fq.gz,fastq.gz}".|
|--sequencing_quality_summary   | Regex pattern for the summary file, default: "sequencing_summary*.txt".|
|--backbone                     | Path to the fasta file containing backbone sequences.|
|--backbone_name                | Name of the sequence to extract from the backbone fasta excluding the starting ">", e.g.:"BB42". |
|--reference                    | Path to the reference genome to use, will ingest all index files in the same directory.|
|--output_dir                   | Directory path where the results, including intermediate files, are stored. |
|--snp_filters.min_dir_ratio, --indel_filters.min_dir_ratio | Minimum ratio of variant-supporting reads in each direction (default: 0.001 (SNP); 0.002 (Indel)).|
|--snp_filters.min_dir_count, --indel_filters.min_dir_count | Minimum number of variant-supporting reads in each direction (default: 5).|
|--snp_filters.min_dpq, --indel_filters.min_dpq             | Minimum positional depth after Q filtering (default: 5_000).|
|--snp_filters.min_dpq_n, --indel_filters.min_dpq_n         | Number of flanking nucleotides to the each position that will determine the window size for local maxima calculation (default = 25).|
|--snp_filters.min_dpq_ratio, --indel_filters.min_dpq_ratio | Ratio of local depth maxima that will determine the minimum depth at each position (default = 0.3).|
|--snp_filters.min_vaf, --indel_filters.min_vaf             | Minimum variant allele frequency (default: 0.003 (SNP); 0.004 (Indel)).|
|--snp_filters.min_rel_ratio, --indel_filters.min_rel_ratio | Minimum relative ratio between forward and reverse variant-supporting reads (default: 0.3 (SNP); 0.4 (Indel)).|
|--snp_filters.min_abq, --indel_filters.min_abq             | Minimum average base quality (default: 70).|


## Using alternative Backbones

Due to lab conditions a different backbone might be used. the --backbone parameter can be set to any fasta file.
The following defaults are available by default in the pipeline, you can enable them by copying the value in the the value column and pasting it behind the cli command. 

|backbone| value |default|
|--------|-------|-------|
|BB22 | --backbone BB22 | |
|BB25 | --backbone BB25 | |
|BB41 | --backbone BB41 | |
|BB42 | --backbone BB42 |X|
|BBCS | --backbone BBCS | |
|BBCR | --backbone BBCR | |



## Roadmap / Todo:
 
 ### Performance:
 1. Multi-threaded variant calling

## changelog

Please see CHANGELOG.md




## Install-dependencies

### Nextflow

Download the latest version by running the example below:

```bash
wget -qO- https://get.nextflow.io | bash
```

or see [The official Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html#installation)

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

for the latest up to date information see their [official documentation](https://apptainer.org/docs/admin/main/installation.html)




## Running on A SLURM cluster such as UMCU HPC

login to the hpc using SSH. there start a sjob with:

```bash
srun --job-name "InteractiveJob" --cpus-per-task 16 --mem=32G --gres=tmpspace:450G --time 24:00:00 --pty bash
```

go to the right project folder

``` bash
cd /hpc/compgen/projects/cyclomics/cycloseq/pipelines/cycloseq/
```

start the pipeline as normal

## Developer notes

### Cycas addition to the repo

Cycas was added as a subtree using code from: https://gist.github.com/SKempin/b7857a6ff6bddb05717cc17a44091202.
This was done instead of submodule to make pulling of the repo easier for endusers and to stay compatible with nextflow run <remote> functionallity.


more specifically:
``` bash
git subtree add --prefix Cycas https://github.com/cyclomics/Cycas 0.4.3 --squash
```

To update run
``` bash
git subtree pull --prefix Cycas https://github.com/cyclomics/Cycas <tag> --squash
```
