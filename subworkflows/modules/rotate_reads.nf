
nextflow.enable.dsl=2


process RotateByCigar{
    // use custom code by liting to covert the reads from their CAB-BCA to ABC
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'


    input:
        tuple path(bam), path(bai)

    output:
        path "*_rotated.fasta"

    script:
        """
        python /rotate_by_cigar.py -i $bam -o ${bam.SimpleName}_rotated.fasta
        """
}
