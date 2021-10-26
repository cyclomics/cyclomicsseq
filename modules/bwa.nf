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
        path fasta
        path ref

    output:
        path "${fasta.simpleName}.bam"
        
    script:
        """
        bwa mem -M -t ${task.cpus} -c 100 $ref $fasta | \
        sambamba view -S -f bam /dev/stdin | \
        sambamba sort -t ${task.cpus} /dev/stdin -o "${fasta.simpleName}.bam"
        """
}
