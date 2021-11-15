
nextflow.enable.dsl=2

process TrimFasta {
    // sed call to remove all info after first comma in every line.
    // used to reduce the length of the readname in fasta files
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    container "ubuntu:20.04"
    
    input:
        each path(fasta)

    output:
        path "*_trimmed.fasta"

    script:
        """
        sed 's/,.*//' $fasta > ${fasta.SimpleName}_trimmed.fasta
        """
}


process ConcatenateFasta {
    // Call `cat` on all files (eg fasta) that enter the process
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    container "ubuntu:20.04"
    
    input:
        path(fasta)

    output:
        path "concat.fasta"

    script:
        """
        cat *.fasta > concat.fasta
        """
}

process Maf2sam {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    container "lastal"

    input:
        path(maf)
    
    output:
        path "${maf.SimpleName}.sam"

    script:
    """
    /lastal_921/scripts/maf-convert sam $maf > ${maf.SimpleName}.sam
    """
}