
nextflow.enable.dsl=2


include {
    FastqToFasta
    MergeFasta
    Extract5PrimeFasta
    Extract3PrimeFasta
    ExtractSpecificRead
} from "./modules/seqkit"

include {
    BwaIndex
    BwaMemReferenceNamedBam as BwaMemReverse
} from "./modules/bwa"

include {
    SamtoolsFlagstatsMapPercentage
    RemoveUnmappedReads
    PrimaryMappedFilter
    SamToBam
    SamtoolsIndexWithID
} from "./modules/samtools"

include {
    FilterBams
} from "./modules/utils"

include{
    Tidehunter53QualTable
    TideHunterFilterTableStartpos
    TideHunterQualTableToFastq
    TideHunterQualTableToJson
    TideHunterQualJsonMerge
} from "./modules/tidehunter"


include {
    Minimap2Index
    MinimapAlignMany
    Minimap2AlignAdaptive
    Minimap2AlignAdaptiveParameterized
} from "./modules/minimap"

include {
    Cycas
    CycasSplit
} from "./modules/cycas"

include {
    MedakaSmolecule
} from "./modules/medaka"

workflow  ConsensusBasic{
    take:
        read_fq_ch
        backbone_fq_ch
    emit:
        Consensus.output
    main:
        Consensus(read_fq_ch)   
}


workflow  ReverseMapping{
    // ReverseMapping aligns all backbones to the read, very IO intensive.
    take:
        read_fastq
        backbone_fasta

    emit:
        Consensus.output
    main:
        FastqToFasta(read_fastq)
        // SplitSequences(fastas)
        BwaIndex(FastqToFasta.out.splitFasta( by: 1 , file: true ))
        BwaMemReverse(backbone_fasta, BwaIndex.out.map { file -> tuple(file.baseName, file)})
        // aligned_reads = BwaMemReverse.out.map { file -> tuple(file.baseName, file)}
        // BwaMemReverse.out.view()
        SamtoolsFlagstatsMapPercentage(BwaMemReverse.out)
        // SamtoolsFlagstatsMapPercentage.out.join(BwaMemReverse.out).view()
        FilterBams(SamtoolsFlagstatsMapPercentage.out.join(BwaMemReverse.out))

}


workflow TidehunterBackBoneQual{
    // TidehunterBackBoneQual takes the backbones and runs tidehunter whilst providing the 3 and 5 prime regions to tidehunter,
    // TODO: rotate?

    take:
        read_fastq
        reference_genome
        backbone_fasta
        backbone_primer_len
        backbone_name
    main:
        // get backbones
        ExtractSpecificRead(backbone_fasta, backbone_name)
        Extract5PrimeFasta(ExtractSpecificRead.out, backbone_primer_len)
        Extract3PrimeFasta(ExtractSpecificRead.out, backbone_primer_len)

        Tidehunter53QualTable(read_fastq.combine(Extract5PrimeFasta.out).combine(Extract3PrimeFasta.out))
        TideHunterFilterTableStartpos(Tidehunter53QualTable.out)
        TideHunterQualTableToFastq(TideHunterFilterTableStartpos.out)
        TideHunterQualTableToJson(TideHunterFilterTableStartpos.out)

        id = TideHunterQualTableToJson.out.first()map( it -> it[0])
        id = id.map(it -> it.split('_')[0])
        jsons = TideHunterQualTableToJson.out.map(it -> it[1]).collect()
        TideHunterQualJsonMerge(id, jsons)
        // BWaMemSorted(Cutadapt.out, reference_gen)
        // SambambaSortSam(BwaMemSorted.out)
        // SamtoolsIndex(SambambaSortSam.out)
        // RotateByCigar(SamtoolsIndex.out)

    emit:
        fastq = TideHunterQualTableToFastq.out
        json = TideHunterQualJsonMerge.out

}

workflow CycasConsensus{
    take:
        read_fastq
        reference_genome
        backbone_fasta
    main:

        Minimap2AlignAdaptiveParameterized(read_fastq, reference_genome)
        SamtoolsIndexWithID(Minimap2AlignAdaptiveParameterized.out)
        Cycas(SamtoolsIndexWithID.out)

    emit:
        fastq = Cycas.out.map( it -> it.take(2))
        id = Cycas.out.first().map( it -> it[0])
        json = id.combine(Cycas.out.map( it -> it[2]))
}

workflow CycasMedaka{
    take:
        read_fastq
        reference_genome
        backbone_fasta
    main:

        MinimapAlignMany(read_fastq, reference_genome)
        SamToBam(MinimapAlignMany.out)
        CycasSplit(SamToBam.out)
        MedakaSmolecule(CycasSplit.out.flatten())

    emit:
        id = CycasSplit.out.first().map( it -> it[0])       
        fastq = MedakaSmolecule.out
        json = id.combine(CycasSplit.out.map( it -> it[2]))

}
