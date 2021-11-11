


process SambambaSortSam{
    publishDir "$baseDir/data/out/$workflow.runName/sambamba/sort"
    container 'mgibio/dna-alignment:1.0.0'
    
    cpus = 6

    input:
        each path(sam)

    output:
        path "*sorted.bam"
        
    script:
        """
        sambamba view -S -f bam $sam | sambamba sort -t ${task.cpus} /dev/stdin -o "${sam.simpleName}_sorted.bam"
        """
}