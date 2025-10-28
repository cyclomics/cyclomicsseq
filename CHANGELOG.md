# changelog

# 1.1.2
- Updated Docker container to version 0.9.1, to include Cycas 0.5.4, which handles an issue where some consensus blocks have invalide quality probabilities. If a read contains invalid quality probabilities, then no consensus is provided for that read.
- Fixed an issue where printing the report fails because no SNPs were found for the scatterplots. If no reads were found, the variant support tab will display a message stating this, and the report HTML will still be produced.


## 1.1.1

- Updated Docker container to version 0.9.0, to include Cycas 0.5.2
- Fixes metadata pairing with consensus BAM files
- Handles spontaneous connection errors when requesting cosmic legacy IDs, now leading to 'None' annotations in case of error

## 1.1.0

- Added variant support scatterplots to the output report.
- Updated Cycas to 0.5.2, allowing the YN tag to be added to the consensus BAM files, encoding consensus read structure.

## 1.0.0

- Major overhaul to CyclomicsSeq, with a simpler workflow architecture. The EPI2ME interface and user guidelines on the `README.md` were also updated.
- Changed naming of important output files to contain the sample ID, a randomly generated unique run ID and the barcode (if present). Users can define the sample ID -- if undefined, a sample ID will be retrieved from a valid parent folder name.
- Added an argument enabling users to provide a BED file with genomic regions to be analysed (`--region_file [regions.bed]`)
- Added support for input directories containing several barcode subdirectories with FASTQ files, where each subdirectory name is interpreted as a sample name. CyclomicsSeq will analyse each sample separately and output separate files and reports per sample.
- Added splitting of large FASTQ files per number of reads (`--split_fastq_by_size [true]` and `--max_fastq_size [40000]`).
- Added a filter for alignment rate to provided or detected regions of interest (`--filter_by_alignment_rate [false]` and `--min_align_rate [0.8]`).
- Updated docker container with seqkit version 2.9.0 and moved container to `cyclomics/cyclomicsseq:0.8.2`.
- Fixed an error in read structure donut plotting where reads with a single alignment were not being counted; adjusted function to be valid with coordinate-sorted BAM files.
- Fixed an issue where alignment was running out of memory due to inappropriate allocation.
- Fixed an issue where the annotated VCF output file was invalid due to a whitespace in COSMIC legacy ID tags.
- Minor plotting and report updates.
- Known contaminant detection has been temporarily disabled due to unknown Freebayes error.

## 0.12.2

- Enable overwriting of nextflow timeline and report
- Hotfix for missing info tab in VCF entries for SNPs in ends of amplicons with some read depth but zero alternative nucleotide counts
- Hotfix for EnsemblVEP API responses that caused a JSON decoding error

## 0.12.1

- Fix an issue where the reference nucleotide couldn't be read due to pysam.FastaFile locking the file, reference nucleotide appeared as '.' instead
- Fix an issue where an alternative allele supported by all reads was not counted, because the total number of observed alleles was not > 1, either overall or in a single read direction

## 0.12.0

- Contaminant detection module
- Make perbase max_depth a CLI option
- Add promethion profile
- Atomize variant filtering, sorting and merging

## 0.11.0

- Added GPLv3 license
- Added a hover tooltip to the alignment count figure in the report HTML
- Fixed calculation of amplicon/region of interest edges, resolving issues with indel detection
- Updated Cycas to version 0.4.8
- Minor refactoring

## 0.10.1

- Fix version number shown in EPI2ME

## 0.10.0

- Changes to EPI2ME interface, improving advanced variant calling options
- Substitute deprecated pd.lookup for pd.factorize while parsing perbase tables
- Changes to HTML report, where plotting error messages are shown as a warning.
- If a report plot fails, the final report is still shown, with a warning message.
- Activate dynamic variant filtering for variant calling potion `validate`
- Ammended how the AO tag is written in vc validate: in accordance with Freebayes, the tag stands for the sum of qualities of the alternative allele
- Streamlining of reporting and variant calling subworkflows

## 0.9.4

- VCF headers have been corrected for validate
- Pipeline version updated on nextflow.config
- Nextflow run reports have been renamed to avoid confusion

## 0.9.3

- Updates to metadata plotting, making it less resource-intensive and including optional subsampling (turned on by default)
- Added Freebayes variant calling with option --metadata.subsample_size flag. if is set to 0, then all reads are taken into account for plotting, essentially inactivating the subsampling feature.
- When using Freebayes variant calling, variant filtering uses a dynamic strategy to filter on minimum VAF (not yet implemented for validate variant calling)
- The output directory for the nxf trace report was adjusted to within the project output directory
- Changes to conda environment and docker container, now including freebayes=1.3.2
- Minor fixes

## 0.9.2

- Changed Region of interest detection
- Improved report, added standard and detailed versions of the report

## 0.9.1

- Changed the mounting of the data in Singularity containers, so that they are compatible with SLURM.
- Changed the locking of files in multithreaded variant calling, to prevent issues in single read file systems.
- Change the allocation of tmp directories in variant calling processing steps.
- Disabled containers for processes that can be executed in the main cyclomicsseq container.
- Disabled the split_on_adapter_sequence by default.
- Changed the default backbone to BBCS, the backbone shipped with the developer access kit.
- Changed the default behavior of overwriting the Nextflow trace file.

## 0.9.0

- Added integration with EPI2ME Labs.
- Changed consensus insert filtering rule to pass reads with YM >= 3 instead of YM > 3.
- Fixed filtered consensus insert number propagation to report. Workflow now outputs BAM with consensus inserts after filtering.
- Fixed VCF position reporting when depth at position is 0
- Fixed non-contiguous chromosome error by adding back -a flag to bcftools sort
- Changed AnntotateBamXTags label
- Added the parameterized roi settings
- Fixed reference nucleotide reporting if majority nucleotide
- Fixed metadata json plotting when no reads where found
- Fixed memory hogging of the donut plot generation by looping over a sorted bam iso making a dataframe
- Added BB41T (legacy development backbone)
- Made the reporting truely optional.
- Fixed bug where very low depth variants could crash the variant filter.
- Fixed many smaller issues.

## 0.8.2

- Fixed issue where variant filtering throwed an error when variant dataframe becomes empty before DPQ filtering.

## 0.8.1

- Default minimum VAF filter for SNPs was lowered to 0.003.
- Docstrings added to variant calling/filtering/reporting functions as needed.
- Documentation on new variant calling user arguments added to README.
- Removed workflows folder.

## 0.8.0

- Major overhaul of variant calling, with new default parameters that are also settable by user.
- Changed Default post consensus alignment strategy to BWA mem.
- FIxed bug where count.txt files where empty.
- minor changes to CI pipeline.

## 0.7.3

- Added correct labels to BWA processes to prevent 137 errors in conda runs.

## 0.7.2

- Added bwa as option in post consensus alignment step
- fixed issues with Conda enviroment in some cases.
- standardized the process labels

## 0.7.1

- Added indel detection
- Added variant effect prediction for Grch38.
- Added visual improvements to the report
- Fixed warning with respect to internal Nextflow tuple
- Moved filtering into a seperate submodule

## 0.7.0

- Change output structure
- Added single page reporting, available in QC/report.html
- Added readsplitting

## 0.6.1

- fix bug in reference genome info generator
- fix conda environment parser
- Add backbone files to the pipeline
- Update README.md

## 0.6.0

- Update cycas version for consensus generation
- Alter variant calling parameters
- Pass vcf parameters to final file
- Changed folder structure

## 0.5.2

- Update Cycas version
- set --min-repeat-count default to 5

## 0.5.0

- Make --minimum-repeat-count parameter settable by user
- resolve tag conflict in both x and y tags
- Add csv output to bam accuracy.
- Increase decimal count in vcf to 6.
- Increased sensitivity of variant calling.

## 0.4.6

- Fix annotation file grouping
- Fix region file code
- Add option to select reads based on repeat count (YM bam tag)

### 0.4.5

- Added Y tags to final bam
- Improved QC

### 0.4.4

- Fix issues with bam accuracy plotting

### 0.4.3

- Fix fontsizes

### 0.4.2

- Added many QC steps with plots,
- Removed MionIONQC
- changed default variant calling to validate
- added auto detect of region of interest

### 0.4.1

- changed default read_pattern regex to detect rebasecalled sequencing runs automatically when pointed at the output folder.

### 0.4.0

- Updated Cycas to prevent runtime error with BB41
- Added variant validation optionality.
- Added quick_results flag for a glance of the results in the terminal.
- Added profiles for conda and singularity support.

### 0.3.1

- Updated cycas version
- disable reporting due to to metadata incompatibility

### 0.3.0

- Added Cycas as the main consensus caller
- Added filtering
- Changed to Varscan variant calling
