# 🌴 Cycas 🌴

CYclomics Consensus Alignment Sequences, or the interesting genus of the ["vredespalm"](https://en.wikipedia.org/wiki/Cycas) EN:(Peace palm).

Cycas uses the alignment against the reference genome to determine the type of read, and creates consensus sequences accordingly.

## TODO

- Better default management. Defaults should only be defined in the "public" space, all the rest should contain no defaults and require explicit.
- Rewrite of the metric system:
  - Rethink of different levels of metrics (per read, per file, per run).
  - Plots per metric, save to a small file ready for plotting. plot/save/load methods per metric.
  - Callback system for metrics to add in app calls.
  - Unit tests and examples for classification/split/metrics/algorithms.
- Rotate consensus reads: complement MSA arrays, such that overlapping nucleotides are identified by some sort of index number.
- Remove deprecated legacy code.

## Changelog

### 0.6.0
 - Major changes to Cycas, enabling splitting of chimeric reads, classification of concatemer structures and calibration of consensus quality scores.
 - Introduced the following commands to `cycas.py`:
   - `consensus` generates consensus from an aligned BAM file, outputs a consensus FASTQ and optionally a metadata JSON.
   - `classify` (for advanced users) splits and classifies reads from an aligned BAM file, outputs a split BAM file, a FASTQ file and optionally an HDF5 file. Does not perform consensus.
   - `simulate` (for advanced users) simulates CyclomicsSeq reads from a reference genome.
   - `file-prep` (for advanced users) prepares data from an HDF5 file for calibration.
   - `fit-model` (for advanced users) fits the model according to the prepared data.
   - `quick-consensus` (DEPRECATED) generates the structure-unaware consensus. Previously `consensus`.
   - `quick-consensus-sr` (DEPRECATED) generates structure-unaware consensus for a single read. Previously `single-read`.
 - Changed segment IDs to be locus-specific, unaware of orientation (if two alignments are in the same region, they will have the same ID irrespective of their orientation).

### 0.5.3
 - Fix an issue in alternative representation where first alignment gap was added twice
 - Fix an issue in alternative representation where last alignment gap was not added
 - Change annotation of negative alignment gap lengths to correctly describe overlaps (e.g. `-6:O` instead of `-6:U`)
 - Add segment start position within concatemer to alternative representation
 - Extend unit tests for alternative representation

### 0.5.2
 - Add YN tags for the consensus structure of reads

### 0.5.1
 - Change user for nextflow compatibility

### 0.5.0
 - Update dependencies to python 3.12 in new docker file
 - Change regex patterns to new syntax
 - Remove unused imports

### 0.4.9
 - Updated docker container
 - Added procps
 - Minor log changes

### 0.4.8
- Optimized `minimal_deletion_ratio` to new optimum (`minimal_deletion_ratio=55`) as described in ELAB experiment 300

### 0.4.7
- Fix Unknown with new SingleInsertUncertain and DoubleInsertUncertain classifications

### 0.4.6
- Fix CI pipeline settings
- Add category for reads with large unmapped repeating segments.
- Add start and stop lines for the read to make read ending less ambiguous.
- Fix spelling mistake in classification.
- Add additional rules for classification

### 0.4.5
- Change distance to detect groups.
- Add requirement that at least 10% of bases need to be present at a location to create consensus.

### 0.4.4
 - fix bug that causes a decrease in the number of insert only reads.
 - Added first CI steps for linting.

### 0.4.3
- Fix in type annotation

### 0.4.2
- Fix issue where 4th gen barcodes cause errors, please note that the barcode reported is incorrect for these in this update.
- Fix bug where extended cigar string was using a preset extender when applying insert locations.

### 0.4.1
- Fix for cases where value cannot be set on missing key.

### 0.4.0
- large overhaul of the codebase in general
- rewrote classifier
- changed metadata structure
- Changed consensus calling code to exclude unstarted and ended sequences.
- fixed bug where cigar strings where not parsed properly in some cases.

### 0.3.0
- Added subcommands `split`, `consensus`, and `single_read`.
- Renamed `Short*` to `LowAlignmentCount*`.
- Changed `TwoInsert` class to .
- Added `Singleton*` classification.
- Added Multimer Backbone classes
- Added `DirectionalFlip` for momomer, renamed old one to `DoubleDirectionalFlip`.
- Added `ComplexConcatemer`
- Pinned colors for each classification for beter comparison.


### 0.2.2
- Added Y tags as defined by cyclomics standards to fastq name line as well as metadata
- Changed warnings with regards to filtering
- Created Barcode extractor
- More general improvements

### 0.2.1

- Added metadata output json option (--create-metadata-json)
- Added classification details output (-create-classification-detail-json)
- Fixed bug with spaced inserts in chromosome

### 0.2.0

- Improved consensus calling
- Added filter at 3 repeats per alignment block
- Tweaked minimum deletion ratio

## Usage

All that is required is an alignment file of the raw Cyclomicseq data, which contains supplementary alignments.
example of usage is:

### Consensus

``` bash
python3 cycas/cycas.py consensus --bam-file <somepath>.bam --outfile <somepath>.fastq --metadata-json <somepath>.json 
```

### Single read

``` bash
python3 cycas/cycas.py single-read --bam-file <somepath>.bam --outfile <somepath>.fastq --read-name <some-read-name>
```

### Options

Options and their availability per subcommand:

|flag                                 |type  |Consensus |Single read |help |
|-------------------------------------|------|----------|------------|-----|
| --bam-file                          | PATH | X        | X          | Input BAM file, indexed. |
| --read-name                         | STR  |          | X          | name of the read as in the BAM file. |
| --output                            | PATH | X        | X          | Filepath for the output file in Fastq format.  |
| --plot-readtypes                    | BOOL | X        | X          | Plots 10 reads, usefulle for testing and classification development.   |
| --metadata-json                     | BOOL | X        | X          | Create a detailed output in json format about the classification  |
| --limit-calls                       | INT  | X        | X          | Limit the number of reads processed, usefull for testing.  |


## Design

First we obtain all alignments, derived from pysam AlignedSegments, called `Alignment` for a read in a `AlignmentGroup`. AlignmentGroups are also called groups.

Cycas uses Rules to determine the class of each read. The hierarchy is :

1. A `Rule`, or `ClassificationRule` is the atomic unit of this codebase. they have a score function that returns values between 0 and 1.
1. There is a `MetadataClassifier`, this group has multiple `Classifier`s it will apply, right now there is only one.
1. The `Classifier` with the lowest priority and that passes all requirements set by the `Rule`s will be reported, this also determines the consensus generation. 

## Classification definitions

1. SingleInsert

    Read where there is only one alignment against the species of origin.

    ![pic of SingleInsert](images/classification_examples/SingleInsert.png)

1. SingleBackBone

    Read where there is only one alignment against a single backbone.

    ![pic of SingleBackBone](images/classification_examples/SingleBackBone.png)

1. BackboneInsert

    A read where we map to the backbone and an insert in an alternating pattern. These reads are as the Cyclomicsseq protocol envisions them.

    ![pic of BackboneInsert](images/classification_examples/BackboneInsert.png)

1. MessyAlignment

    Any read with more than 3 different directional chromosome alignments eg chr17+, chr17- BB22+ and chr2- 

    ![pic of ComplexConcatemer](images/classification_examples/ComplexConcatemer.png)

1. Unkown

    Read that we where not able to classify.

    ![pic of Unkown](images/classification_examples/Unkown.png)

## Colors

The main idea behind the piechart colors is that Blue is good, and orange is bad. The darker the hue of blue the better the read is. The darker the hue of orange, the worse the read is.

### old structure

35:U,
98:I:F:TP53:15160:97:0,
1:U,
248:BB:R:BB25:0:247:0,
134:I:F:TP53:15108:139:0,
5:U

bases_on_read:type:orient:assembly:position:mapping_length:ID
 
