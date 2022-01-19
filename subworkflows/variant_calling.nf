include {
    Freebayes
} from "../modules/freebayes"

include {
    CreateRererenceDict
    GatkAddOrReplaceReadGroups
    Mutect2TumorOnly
} from "../modules/gatk"


workflow FreebayesSimple{
    take:
        reads
        reference_genome

    main:
        Freebayes(reads.combine(reference_genome))
    emit:
        Freebayes.out
}

workflow Mutect2{
    take:
        reads
        reference_genome
    main:
        CreateRererenceDict(reference_genome)
        GatkAddOrReplaceReadGroups(reads)
        Mutect2TumorOnly(GatkAddOrReplaceReadGroups.out.combine(reference_genome).combine(CreateRererenceDict.out))
    emit:
        Mutect2TumorOnly.out
}
