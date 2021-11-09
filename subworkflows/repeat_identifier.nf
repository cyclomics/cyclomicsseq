
nextflow.enable.dsl=2

include{
    DummyProcess as Rotate
} from "../modules/dummy.nf"

// include {
//     Rotate
// } from "../modules/"

include {
    BwaMemSorted
} from "../modules/bwa.nf"

include {
    RotateByCigar
} from "../modules/rotate_reads.nf"

include {
    SamtoolsIndex
} from "../modules/samtools.nf"

include {
    Tidehunter
} from "../modules/tidehunter.nf"

include {
    TrimFasta
    ConcatenateFasta
} from "../modules/utils.nf"

workflow  TideHunterBasic{
    take:
        read_fq_ch
        reference_genome
    emit:
        RotateByCigar
    main:
        Tidehunter(read_fq_ch)
        // collect (process in one go) all the tidehunter files
        ConcatenateFasta(Tidehunter.out.collect())
        TrimFasta(ConcatenateFasta.out)
        BwaMemSorted(TrimFasta.out, reference_genome)
        SamtoolsIndex(BwaMemSorted.out)
        // we need indexed bams to rotate them
        RotateByCigar(SamtoolsIndex.out)
}
