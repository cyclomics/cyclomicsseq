include {
    Freebayes
} from "./modules/freebayes"

include {
    CreateRererenceDict
    GatkAddOrReplaceReadGroups
    Mutect2TumorOnly
} from "./modules/gatk"

include {
    VarscanFiltered
} from "./modules/varscan"

include {
    VcfToBed
} from "./modules/bedops"

include {
    VariantValidate
    FilterVariants
} from "./modules/bin"

include {
    FindRegionOfInterest
} from "./modules/samtools"

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

workflow Varscan{
    take:
        reads
        reference_genome
    main:
        VarscanFiltered(reads.combine(reference_genome))
    emit:
        VarscanFiltered.out
}
workflow ValidatePosibleVariantLocations{
    // Allow to determine VAF for given genomic positions in both bed and vcf format
    take:
        reads_aligned
        variant_file_name
        reference

    main:
        variant_file = Channel.fromPath(params.region_file, checkIfExists: false)

        if (variant_file_name == 'auto') {
            positions = FindRegionOfInterest(reads_aligned)
        }
        else if( variant_file_name.endsWith('.vcf') ) {
            positions = VcfToBed(variant_file)
        }
        else {
            positions = variant_file
        }
        VariantValidate(reads_aligned, positions)
        FilterVariants(VariantValidate.out)

    emit:
        locations = VariantValidate.out
        variants = FilterVariants.out
}
