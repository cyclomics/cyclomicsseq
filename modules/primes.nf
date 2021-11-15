
nextflow.enable.dsl=2


// Steps to extract prime end of fastq files
process Extract5PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    container "staphb/seqtk:1.3"
    
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

process Extract3PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "staphb/seqtk:1.3"

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

