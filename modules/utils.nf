
nextflow.enable.dsl=2

process TrimFasta {
    // sed call to remove all info after first comma in every line.
    // used to reduce the length of the readname in fasta files
    publishDir "$baseDir/data/out/$workflow.runName/trim"

    container "ubuntu:20.04"
    
    input:
        each path(fasta)

    output:
        path "trimmed.fasta"

    script:
        """
        sed 's/,.*//' $fasta > trimmed.fasta
        """
}


process ConcatenateFasta {
    // Call `cat` on all files (eg fasta) that enter the process
    publishDir "$baseDir/data/out/$workflow.runName/concat"

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