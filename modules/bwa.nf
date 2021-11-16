
nextflow.enable.dsl=2

process BwaIndex{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'
    
    input:
        each path(fasta)
        tuple val(sampleId), file(reference)

    output:
        path "${fasta.simpleName}.bam"
        
    script:
        // TODO: add a grep to find out if ref is fa or fasta
        // piplefail for better control over failures
        """
        set -euxo pipefail
        bwa mem -M -t ${task.cpus} -c 100 $reference $fasta | \
        sambamba view -S -f bam /dev/stdin | \
        sambamba sort -t ${task.cpus} /dev/stdin -o "${fasta.simpleName}.bam"
        """
}

process BwaMem16c{
    // Run bwa with sam output using 16 cores
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'
    
    cpus = 16

    input:
        each path(fasta)
        tuple val(sampleId), file(reference)

    output:
        path "${fasta.simpleName}.sam"
        
    script:
        // TODO: add a grep to find out if ref is fa or fasta
        """
        set -euxo pipefail
        bwa mem -M -t ${task.cpus} -c 100 $reference $fasta > ${fasta.simpleName}.sam
        """
}
