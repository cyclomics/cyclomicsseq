# CycloSeq

Pipelines to process Cyclomics data.

This pipeline uses prior information from the backbone to increase the effectiveness of consensus calling from the circular DNA protocol by Cyclomics.  

## Dependencies

Click for installation instructions:

- Nextflow[#install-dependencies]
- Docker[#install-dependencies]
- Acces to the Github repo and a valid PAT token[https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token]


### Data requirements

- data output by ONT Guppy-hac
- Backbone fasta
- Reference genome


## Usage

In this section we assume that you have docker and nextflow installed on your system, if so running the pipeline is straightforward. You can run the pipeline directly from this repo, or pull it yourself and point nextflow towards it.

If you want to run the pipeline directly from github you need to use a Personal Access Token (PAT) as the password. Click the link to see how to create a [Personal Access Token (PAT)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). Your PAT should have at least read permissions.
an example whould be:

```bash
nextflow run cyclomics/cycloseq -user <github_username> -r <current_version> -profile docker --input_read_dir '/sequencing/20220209_1609_X3_FAS06478_0ed4361c' --output_dir '/data/myresults' --reference '/data/reference/chm13v2.fasta' --backbone-fasta '/data/reference/cyclomics_backbones.fasta' 
```

Fill in your PAT when promted for it.



### Singularity

If docker is not an option, singularity(or Apptainer as it is called since Q2 2022) is a good alternative that does not require root access and therefor used in shared compute environments.


### Conda

IF you pull the repo from github for local excecution and want to use conda as the runtime environment you need to provide an extra additional argument:

```bash
--user_conda_location <path to environemnt.yml>
```

this file is located in the repo, the pipeline needs to know where this file is to run with the correct versions of the required software.

This means the full command becomes:
```bash
nextflow run cyclomics/cycloseq -profile conda --user_conda_location /path/to/environemnt.yml --input_read_dir '/sequencing/20220209_1609_X3_FAS06478_0ed4361c' --output_dir '/data/myresults' --reference '/data/reference/chm13v2.fasta' --backbone-fasta '/data/reference/cyclomics_backbones.fasta' 
```



## Roadmap / Todo:
 
 ### Usability:
 1. Implement warning when the backbone is not present in the data
 1. give warning when user sets backbone and its not used




### Flag descriptions

|flag                           | info  |
|-------------------------------|---|
|--input_read_dir               | Directory where the output of Guppy is located, e.g.: "/data/guppy/exp001/".|
|--read_pattern                 | Regex pattern to look for fastq's in the read directory, defaults to: "fastq_pass/**.{fq,fastq}".|
|--sequencing_quality_summary   | Regex pattern for the summary file, default: "sequencing_summary*.txt".|
|--backbone_fasta               | Path to the fasta file containing backbone sequences.|
|--backbone_name                | Name of the sequence to extract from the backbone fasta excluding the starting ">", e.g.:"BB22". |
|--reference | Path to the reference genome to use, will ingest all index files in the same directory.|
|--output_dir | Directory path where the results, including intermediate files, are stored. |

## changelog

Please see CHANGELOG.md

### Running on A SLURM cluster such as UMCU HPC

login to the hpc using SSH. there start a sjob with:

```bash
srun --job-name "InteractiveJob" --cpus-per-task 16 --mem=32G --gres=tmpspace:450G --time 24:00:00 --pty bash
```

go to the right project folder

``` bash
cd /hpc/compgen/projects/cyclomics/cycloseq/pipelines/cycloseq/
```

start the pipeline as normal


## Install-dependencies

### Nextflow

wget -qO- https://get.nextflow.io | bash

### Conda

Download the latest conda version from here https://docs.conda.io/en/latest/miniconda.html#linux-installers

Run the below command and follow process:
```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

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
