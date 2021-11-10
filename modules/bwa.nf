
nextflow.enable.dsl=2

process BwaIndex{
    publishDir "$baseDir/data/out/$workflow.runName/bwa/index"
    container 'biocontainers/bwa:v0.7.17_cv1'
    
    input:
        path reference_genome

    output:
        path "${reference_genome}.*", emit: bwa_index

    script:
        """
        bwa index $reference_genome
        """
}

process BwaMemSorted{
    publishDir "$baseDir/data/out/$workflow.runName/bwa/sorted"
    container 'mgibio/dna-alignment:1.0.0'
    
    input:
        each path(fasta)
        tuple val(sampleId), file(reference)

    output:
        path "${fasta.simpleName}.bam"
        
    script:
        // TODO: add a grep to find out if ref is fa or fasta
        """
        bwa mem -M -t ${task.cpus} -c 100 ${sampleId}.fa $fasta | \
        sambamba view -S -f bam /dev/stdin | \
        sambamba sort -t ${task.cpus} /dev/stdin -o "${fasta.simpleName}.bam"
        """
}
