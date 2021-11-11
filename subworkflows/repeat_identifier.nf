
nextflow.enable.dsl=2

include {
    BwaMemSorted
    BwaMem16c
} from "../modules/bwa.nf"

include {
    RotateByCigar
} from "../modules/rotate_reads.nf"

include {
    Extract5PrimeFasta
    Extract3PrimeFasta
} from "../modules/primes.nf"

include {
    SamtoolsIndex
} from "../modules/samtools.nf"

include {
    SambambaSortSam
} from "../modules/sambamba.nf"

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
    main:
        Tidehunter(read_fq_ch)
        // collect (process in one go) all the tidehunter files
        ConcatenateFasta(Tidehunter.out.collect())
        TrimFasta(ConcatenateFasta.out)
        BwaMem16c(TrimFasta.out, reference_genome)
        SambambaSortSam(BwaMem16c.out)
        SamtoolsIndex(SambambaSortSam.out)
        // we need indexed bams to rotate them
        RotateByCigar(SamtoolsIndex.out)
    emit:
        RotateByCigar.out
}

// workflow TidehunterBackBone {
//     take:
//         read_fq_ch
//         reference_genome
//         backbone
//     main:
//         // get backbones
//         p5_babo = Extract5PrimeFasta(backbone, 10)
//         p3_babo = Extract3PrimeFasta(backbone, 10)

//         Tidehunter53(read_fq_ch, p3_babo, p5_babo)
//         ConcatenateFasta(Tidehunter53.out.collect())
//         TrimFasta(ConcatenateFasta.out)
//         Cutadapt(TrimFasta.out)
//         BWaMemSorted(Cutadapt.out, reference_gen)
//         SambambaSortSam(BwaMemSorted.out)
//         SamtoolsIndex(SambambaSortSam.out)
//         RotateByCigar(SamtoolsIndex.out)

//     emit:
//      RotateByCigar.out
// }
