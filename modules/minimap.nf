

process MinimapAlign{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'staphb/minimap2:2.23'
    cpus = 2
    memory '16GB'

    input:
        tuple val(X), path(fasta), path(reference_genome)

    output:
        tuple val(X), path("${X}.sam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fasta > ${X}.sam
        """
}