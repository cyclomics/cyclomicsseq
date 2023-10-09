include {
    Freebayes
    FilterFreebayesVariants
    SeparateMultiallelicVariants
} from "./modules/freebayes"

include {
    CreateRererenceDict
    GatkAddOrReplaceReadGroups
    Mutect2TumorOnly
} from "./modules/gatk"

include {
    VcfToBed
} from "./modules/bedops"

include {
    FindVariants
    FilterValidateVariants
    MergeNoisyVCF
    MergeFilteredVCF
    AnnotateVCF
} from "./modules/bin"

include {
    FindRegionOfInterest
} from "./modules/samtools"

include {
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
} from "./modules/perbase"

workflow ProcessTargetRegions{
    take:
        variant_file_name
        reads_aligned
    
    main:
        variant_file = Channel.fromPath(params.region_file, checkIfExists: false)

        if (variant_file_name == 'auto') {
            positions = FindRegionOfInterest(reads_aligned)
        }
        else if( variant_file_name.endsWith('.bed') ) {
            positions = variant_file
        }
        else if( variant_file_name.endsWith('.vcf') ) {
            positions = VcfToBed(variant_file)
        }
        else {
            positions = variant_file.map(it -> tuple(it.SimpleName, it))
        }
    
    emit:
        positions

}

workflow CallVariantsFreebayes{
    take:
        reads_aligned
        positions
        reference

    main:
        Freebayes(reads_aligned.combine(reference), positions)
        SeparateMultiallelicVariants(Freebayes.out)
        PerbaseBaseDepthConsensus(reads_aligned.combine(reference), positions, 'consensus.tsv')
        FilterFreebayesVariants(SeparateMultiallelicVariants.out.combine(PerbaseBaseDepthConsensus.out))
        AnnotateVCF(FilterFreebayesVariants.out)

    emit:
        locations = Freebayes.out
        variants = AnnotateVCF.out
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


workflow ValidatePosibleVariantLocations{
    // Allow to determine VAF for given genomic positions in both bed and vcf format
    take:
        reads_aligned
        positions
        reference

    main:
        FindVariants(reference, reads_aligned, positions)
        PerbaseBaseDepthConsensus(reads_aligned.combine(reference), positions, 'consensus.tsv')
        FilterValidateVariants(FindVariants.out.combine(PerbaseBaseDepthConsensus.out))
        MergeNoisyVCF(FindVariants.out)
        MergeFilteredVCF(FilterValidateVariants.out)
        AnnotateVCF(MergeFilteredVCF.out)

    emit:
        locations = MergeNoisyVCF.out
        variants = AnnotateVCF.out
}
