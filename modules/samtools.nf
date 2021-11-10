
nextflow.enable.dsl=2


process SamtoolsIndex{
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        each path(bam)

    output:
        tuple path(bam), path("*.bai")

    script:
        """
        samtools index $bam
        """

}

