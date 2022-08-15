
process FindRegionOfInterest{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path("${bam.simpleName}.thresholds.bed.gz"), path("${bam.simpleName}.roi.tsv")
    
    script:
        """
        mosdepth -n --threads 10 --by 1000 --thresholds 10,100,500,1000,5000 ${bam.simpleName} $bam
        zcat ${bam.simpleName}.thresholds.bed.gz | awk '\$5 > 0' > ${bam.simpleName}.roi.tsv 
        """
}

