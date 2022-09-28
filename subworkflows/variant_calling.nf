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
    FindSNPs
    FilterSNPs
    FindIndels
    MergeVCF
    AnnotateVCF
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
            
            positions = variant_file.map(it -> tuple(it.SimpleName, it))

        }
        FindSNPs(reads_aligned, positions)
        FilterSNPs(FindSNPs.out)
        FindIndels(reference, reads_aligned, positions)
        MergeVCF(FilterSNPs.out.combine(FindIndels.out))
        AnnotateVCF(MergeVCF.out)


    emit:
        locations = VariantValidate.out
        variants = AnnotateVCF.out
}
