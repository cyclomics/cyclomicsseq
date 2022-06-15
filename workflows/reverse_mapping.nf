#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/Reverse Calling
========================================================================================
    Github : tbd
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.input_fastq_files = "/media/dami/cyclomics_003/raw_data/ENZA_0003/raw_data/MAR6933/fastq_pass/FAR67932_pass_f59d80e7_[258,259].fastq.gz"
params.target = "/media/dami/cyclomics_003/raw_data/ENZA_0003/barcodes/ENZA_IDs.fa"
params.output_dir = "data/$workflow.runName"

// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CycloSeq : Cyclomics consensus pipeline
    ===================================================
        input_fastq_files : $params.input_fastq_files
        target            :$params.target
        output_dir        :$params.output_dir
"""


process FastqToFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    cpus 1

    input:
        path(fastq)
    output:
        path("${fastq.SimpleName}.fa")

    script:
        """
        seqkit fq2fa $fastq -o ${fastq.SimpleName}.fa
        """
}

process splitSequences {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(fasta)

    output:
        path("${fasta.SimpleName}_seq_*.fa")

    script:
        """
        awk '/^>/{f="${fasta.SimpleName}_seq_"++d".fa"} {print > f}' < $fasta
        """
}

process BwaIndex{
    // dont publish as coping takes to long
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/bwa:v0.7.17_cv1'
    input:
        path(reference_genome)

    output:
        path"${reference_genome}*", includeInputs: true, emit: bwa_index

    script:
        """
        bwa index $reference_genome
        """
}
process BwaMemReverse{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'

    cpus = 1 

    input:
        each path(fasta)
        path(reference)

    output:
        path "*.sam"
        
    script:
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        """
        REF_START=\$(ls | grep -E ".ann" | cut -d. -f1)
        REF=\$(ls | grep -E "\$REF_START.(fasta\$|fa\$)")
        bwa mem -M -t ${task.cpus} \$REF $fasta > \$REF_START.sam
        """
}

process SamtoolsFilterMapped{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    cpus 1

    input:
        path(bam)

    output:
        path("*.map.bam") optional true

    script:
        """
        for ((i=46; i>=35; i--))
            do
                echo \$i
                samtools view -b -e '[AS] >'\$i $bam > ${bam.SimpleName}_\$i.map.bam
                if (( \$(samtools view ${bam.SimpleName}_\$i.map.bam | wc -l) > 0 )); then
                    break
                fi
                rm ${bam.SimpleName}_\$i.map.bam
            done

        """
}

workflow{
    fastqs = Channel.fromPath(params.input_fastq_files, checkIfExists: true)
    target = Channel.fromPath(params.target, checkIfExists: true)

    fastas = FastqToFasta(fastqs)
    fastas_split = splitSequences(fastas)
    indexed_fasta_reads = BwaIndex(fastas_split.flatten())
    reverse_mapped_reads = BwaMemReverse(target, indexed_fasta_reads)
    filtered_mapped_reads = SamtoolsFilterMapped(reverse_mapped_reads)
}