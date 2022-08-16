# changelog

## 0.6.1
- fix bug in reference genome info generator
- fix conda environment parser

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

