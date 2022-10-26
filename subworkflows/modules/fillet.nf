

process SplitReadsOnAdapterSequence {
    // https://github.com/nanoporetech/duplex-tools/blob/master/fillet.md

     label 'many_cpu_intensive'


    input:
        path(fastq)

    output:
        tuple path("results/${fastq.SimpleName}_split.fastq.gz")

    script:
        """
        duplex_tools split_on_adapter . results/ Native 
        """
} 
