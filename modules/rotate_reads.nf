
nextflow.enable.dsl=2


process RotateByCigar{
    // use custom code by liting to covert the reads from their CAB-BCA to ABC
    publishDir "$baseDir/data/out/$workflow.runName/rotate/"
    container 'cigar_rotate:0.0.1'

    input:
        tuple path(bam), path(bai)

    output:
        path "*_rotated.fasta"

    script:
        """
        python /rotate_by_cigar.py -i $bam -o ${bam.SimpleName}_rotated.fasta
        """
}
