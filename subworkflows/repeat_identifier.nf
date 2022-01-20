
nextflow.enable.dsl=2

include {
    BwaMemSorted
    BwaMem16c
    BwaMem
    BwaMem as BwaMemRealign
} from "./modules/bwa"

include {
    RotateByCigar
} from "./modules/rotate_reads"

include {
    Extract5PrimeFasta
    Extract3PrimeFasta
} from "./modules/seqkit"

include {
    SamtoolsIndex
    SamtoolsSort
    SamtoolsMerge
} from "./modules/samtools"

include {
    SambambaSortSam
} from "./modules/sambamba"

include {
    Tidehunter
    TidehunterLongest
    TideHunterTrimmmer
    TideHunterTableToFasta
} from "./modules/tidehunter"

include {
    TrimFasta
    ConcatenateFasta
} from "./modules/utils"


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
        consensus_tsv = TidehunterLongest(read_fq_ch)
        consensus_fasta = TideHunterTableToFasta(consensus_tsv)
        mapped_consensus = BwaMem(consensus_fasta, reference_genome)

        sorted_mapped_consensus = SamtoolsSort(mapped_consensus)

        // TODO make an append like rolling process
        single_bam = SamtoolsMerge(sorted_mapped_consensus.collect(), "abc")

        // TODO debug rotateby cigar

        SamtoolsIndex(single_bam)
        RotateByCigar(SamtoolsIndex.out)
        // BwaMemRealign(RotateByCigar.out, reference_genome)
    emit:
        RotateByCigar.out
}


workflow TideHunterQuality{
    take:
        read_fq_ch
        reference_genome

}

workflow TidehunterBackBone {
    take:
        read_fq_ch
        reference_genome
        backbone
    main:
        // get backbones
        // TODO: get backbone from a fasta so no custom work is needed
        // TODO: or pick the most occuring backbone 
        backbone_fa = ExtractRead(backbone, "BB22")
        p5_babo = Extract5PrimeFasta(backbone, params.tidehunter.primer_length)
        p3_babo = Extract3PrimeFasta(backbone, params.tidehunter.primer_length)

        Tidehunter53(read_fq_ch, p3_babo, p5_babo)
        ConcatenateFasta(Tidehunter53.out.collect())
        TrimFasta(ConcatenateFasta.out)
        Cutadapt(TrimFasta.out)
        BWaMemSorted(Cutadapt.out, reference_gen)
        SambambaSortSam(BwaMemSorted.out)
        SamtoolsIndex(SambambaSortSam.out)
        RotateByCigar(SamtoolsIndex.out)

    emit:
     RotateByCigar.out
}
