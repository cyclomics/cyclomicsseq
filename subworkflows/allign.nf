
nextflow.enable.dsl=2


include {
    BwaMemSorted
} from "../modules/bwa.nf"


workflow AlignBWA{
    take:
        read_fq_ch
        reference_genome
    emit:
        BwaMemSorted
    main:
        BwaMemSorted(read_fq_ch, reference_genome)
}
