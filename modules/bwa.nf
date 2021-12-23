
nextflow.enable.dsl=2

process BwaIndex{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/bwa:v0.7.17_cv1'
    
    input:
        path reference_genome

    output:
        path "${reference_genome}*", includeInputs: true, emit: bwa_index

    script:
        """
        bwa index $reference_genome
        """
}

process BwaMem{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'

    cpus = 2

    input:
        each path(fasta)
        tuple val(sampleId), file(reference)

    output:
        path "${fasta.simpleName}.sam"
        
    script:
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        """
        REF=\$(ls | grep -E "${sampleId}.(fasta\$|fa\$)")
        bwa mem -M -t ${task.cpus} -c 100 \$REF $fasta > ${fasta.simpleName}.sam
        """
}

process BwaMemReferenceNamedBam{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'

    cpus = 1

    input:
        each path(fasta)
        tuple val(sampleId), file(reference)

    output:
        tuple val("${sampleId.first()}"), path("${sampleId.first()}.bam") , path("${sampleId.first()}.bam.bai")
        
    script:
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        """
        REF=\$(ls | grep -E "${sampleId}.(fasta\$|fa\$)")
        bwa mem -t ${task.cpus} -c 100 -M  \$REF $fasta | sambamba view -S -f bam /dev/stdin | sambamba sort -o ${sampleId.first()}.bam /dev/stdin
        sambamba index ${sampleId.first()}.bam
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
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        // piplefail for better control over failures
        """
        REF=\$(ls | grep -E "${sampleId}.(fasta\$|fa\$)")
        set -euxo pipefail
        bwa mem -M -t ${task.cpus} -c 100 \$REF $fasta | \
        sambamba view -S -f bam /dev/stdin | \
        sambamba sort -t ${task.cpus} /dev/stdin -o "${fasta.simpleName}.bam"
        """
}

process BwaMem16c{
    // Run bwa with sam output using 16 cores
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'
    
    cpus = 3

    input:
        each path(fasta)
        tuple val(sampleId), file(reference)

    output:
        path "${fasta.simpleName}.sam"
        
    script:
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        """
        REF=\$(ls | grep -E "${sampleId}.(fasta\$|fa\$)")
        bwa mem -M -t ${task.cpus} -c 100 \$REF $fasta > ${fasta.simpleName}.sam
        """
}
