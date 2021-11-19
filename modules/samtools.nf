
nextflow.enable.dsl=2


process SamtoolsIndex{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
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

process SamtoolsSort{
    // Given a sam or bam file, make a sorted bam
    // Does break when sam is not `propper` eg no @SQ tags
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        path(input_sam)

    output:
        path "${input_sam.SimpleName}_sorted.bam"

    script:
        """
        samtools sort $input_sam -O bam -o "${input_sam.SimpleName}_sorted.bam"
        """

}

process SamtoolsMerge{
    //  samtools merge – merges multiple sorted files into a single file 

    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        path(input_sam)
        val output_base

    output:
        path "${output_base}.bam"

    script:
        """
        samtools merge -O bam ${output_base}.bam $input_sam
        """

}


process SamtoolsQuickcheck{
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        path(input_sam)

    output:
        stdout

    script:
        """
        samtools quickcheck $input_sam && echo 'Samtools quickcheck ok' || echo 'Samtools quickcheck fail!'
        """
}

process SamtoolsFlagstats{
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        path(input_sam)
    
    output:
        stdout

    script:
        """
        samtools flagstat $input_sam
        """
}
