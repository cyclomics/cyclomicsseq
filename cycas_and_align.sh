# exit when any command fails
set -e

# use a timestamp for tracking
test_id=$(date +%s)

InputBam='~/data/LSK109_26/reads.bam'
outfastq="~/data/LSK109_26/reads_${test_id}.fastq"
outbam="~/data/LSK109_26/reads_${test_id}.bam"
reference='~/data/references/Homo_sapiens/GRCh38/backbone/GRCh38_full_analysis_set_plus_decoy_hla.fa'

# CycasOptions=' --create-classification-detail-json --plot-piechart --limit-calls 50'
CycasOptions=' --create-classification-detail-json --plot-piechart'

# echo $test_id
# echo $InputBam
# echo $outfastq
# echo $outbam
# echo $reference
# echo $CycasOptions

python3 cycas/cycas.py --bam-file $InputBam --output $outfastq $CycasOptions
minimap2 -ax map-ont -t 8 $reference $outfastq  | samtools sort -o $outbam -
samtools index $outbam

echo $test_id
echo $outbam

