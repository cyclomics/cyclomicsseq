#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FastqToFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'pegi3s/seqkit:2.1.0'
    label 'many_cpu_medium'

    input:
        path(fastq)

    output:
        path("${fastq.SimpleName}.fa")

    script:
        """
        seqkit fq2fa $fastq -o ${fastq.SimpleName}.fa
        """
}

process Extract5PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "staphb/seqtk:1.3"
    label 'many_cpu_medium'

    input:
        path fasta
        val length

    output:
        path "*_5p.fasta"

    script:
        """
        seqtk trimfq -L $length $fasta > ${fasta.simpleName}_${length}_5p.fasta
        """
}

process MergeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "staphb/seqtk:1.3"
    label 'many_cpu_medium'

    input:
        path fasta1
        path fasta2
    
    output:
        path "${fasta1.simpleName}_${fasta2.simpleName}.fasta"
    
    script:
        """
        cat $fasta1 > ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        cat $fasta2 >> ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        """

}

process Extract3PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "staphb/seqtk:1.3"
    label 'many_cpu_medium'

    input:
        path fasta
        val length

    output:
        path "*_3p.fasta"

    script:
    // flip, trim, flip back
        """
        seqtk seq -r $fasta | seqtk trimfq -L $length - | seqtk seq -r - > ${fasta.simpleName}_${length}_3p.fasta
        """
}

process ExtractSpecificRead{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'pegi3s/seqkit:2.1.0'
    label 'many_cpu_medium'

    input:
        path fasta
        val readname

    output:
        path "*_${readname}.fasta"

    script:
        """
        seqkit grep -p $readname  $fasta  > ${fasta.simpleName}_${readname}.fasta
        """
}
