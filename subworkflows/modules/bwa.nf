#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BwaIndex{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

    input:
        path reference_genome

    output:
        path "${reference_genome}*", emit: bwa_index

    script:
        """
        bwa index $reference_genome
        """
}

process BwaMem{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

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
        bwa mem -M -t ${task.cpus} -c 100 -L  \$REF $fasta > ${fasta.simpleName}.sam
        """
}

process BwaMemReferenceNamedBam{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

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
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

    input:
        each path(fastq)
        file(reference)
        file(reference_indexes)

    output:
        tuple val("${fastq.simpleName}"), path("${fastq.simpleName}.bam"), path("${fastq.simpleName}.bam.bai")
        
    script:
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        // piplefail for better control over failures
        // REMOVED. TODO: check reference file itself
        """
        bwa mem -R "@RG\\tID:${params.bwamem.readgroup}\\tSM:${params.bwamem.sampletag}\\tPL:${params.bwamem.platform}" -M -t ${task.cpus} -c ${params.bwamem.mem_max_genome_occurance} -L ${params.bwamem.softclip_penalty} -M $reference $fastq | \
        samtools sort -@ ${task.cpus} /dev/stdin -o "${fastq.simpleName}.bam"
        samtools index ${fastq.simpleName}.bam
        """
}

process BwaMemContaminants{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_med_cpu_huge_mem'

    input:
        each path(fastq)
        file(reference)
        file(reference_indexes)

    output:
        tuple val("${fastq.simpleName}"), path("${fastq.simpleName}_contaminants.bam"), path("${fastq.simpleName}_contaminants.bam.bai")
        
    script:
        """
        bwa mem -R "@RG\\tID:${params.bwamem.readgroup}\\tSM:${params.bwamem.sampletag}\\tPL:${params.bwamem.platform}" -M -t ${task.cpus} -c ${params.bwamem.mem_max_genome_occurance} -L ${params.bwamem.softclip_penalty} -M $reference $fastq | \
        samtools sort -@ ${task.cpus} /dev/stdin -o "${fastq.simpleName}_contaminants.bam"
        samtools index ${fastq.simpleName}_contaminants.bam
        """
}

process BwaMem16c{
    // Run bwa with sam output using 16 cores
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // Legacy process
    // container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

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
