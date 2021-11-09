
nextflow.enable.dsl=2


include {DummyProcess as RepeatSplit;
} from "../modules/dummy.nf"

include {
    
} from "../modules/bwa.nf"


workflow RepeatSplitBasic{
    take:
        read_fq_ch
        reference_genome
    emit:
        BwaMemSorted
    main:
        BwaMemSorted(read_fq_ch, reference_genome)
}
