
nextflow.enable.dsl=2

include {
    BwaMemSorted
    BwaMem16c
    BwaMem
} from "../modules/bwa"

include {
    RotateByCigar
} from "../modules/rotate_reads"

include {
    Extract5PrimeFasta
    Extract3PrimeFasta
} from "../modules/primes"

include {
    SamtoolsIndex
    SamtoolsSort
    SamtoolsMerge
} from "../modules/samtools"

include {
    SambambaSortSam
} from "../modules/sambamba"

include {
    Tidehunter
    TidehunterLongest
    TideHunterTrimmmer
} from "../modules/tidehunter"

include {
    TrimFasta
    ConcatenateFasta
} from "../modules/utils"


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

workflow TideHunterKeepLongest{
    take:
        read_fq_ch
        reference_genome
    main:
        consensus = TidehunterLongest(read_fq_ch)
        trimmed_consensus = TideHunterTrimmmer(consensus)
        // TODO make an append like rolling process
        mapped_consensus = BwaMem(trimmed_consensus.trimmed_fasta, reference_genome)
        sorted_mapped_consensus = SamtoolsSort(mapped_consensus)
        
        single_bam = SamtoolsMerge(sorted_mapped_consensus.collect(), "abc")
        // TODO debug rotateby cigar

        // SamtoolsIndex(single_bam)
        // RotateByCigar(SamtoolsIndex.out)
    emit:
        single_bam
}


workflow TideHunterQuality{
    take:
        read_fq_ch
        reference_genome

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
