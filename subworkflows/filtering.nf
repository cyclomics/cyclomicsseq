
nextflow.enable.dsl=2

include {
    FilterShortReads
} from "./modules/seqkit"

include {
    SplitReadsOnAdapterSequence
} from "./modules/fillet"

workflow FilterWithAdapterDetection{
    take:
        read_fq_ch

    main:
        if (params.split_on_adapter != "no") {
            fastq = SplitReadsOnAdapterSequence(read_fq_ch)
        }
        else {
            fastq = read_fq_ch
        }
        FilterShortReads(fastq)
        
    emit:
        FilterShortReads.out
}