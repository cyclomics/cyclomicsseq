

process SplitReadsOnAdapterSequence {
    // https://github.com/nanoporetech/duplex-tools/blob/master/fillet.md
    label 'many_low_cpu_low_mem'

    input:
        path(fastq)

    output:
        path("results/${fastq.SimpleName}_split.fastq.gz")

    script:
        """
        duplex_tools split_on_adapter . results/ Native 
        """
} 
