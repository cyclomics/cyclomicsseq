include {
    Freebayes
    FreebayesContaminants
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
    SortVCF as SortNoisySnp
    SortVCF as SortNoisyIndel
    SortVCF as SortFilteredSnp
    SortVCF as SortFilteredIndel
    MergeVCF as MergeNoisyVCF
    MergeVCF as MergeFilteredVCF
    IntersectVCF
    FilterValidateVariants as FilterSnp
    FilterValidateVariants as FilterIndel
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

workflow CallContaminantMutants{
    take:
        reads_aligned
        positions
        reference

    main:
        FreebayesContaminants(reads_aligned.combine(reference), positions)
        SeparateMultiallelicVariants(FreebayesContaminants.out)
        // AnnotateVCF(FilterFreebayesVariants.out)

    emit:
        locations = FreebayesContaminants.out.map(it -> it[0, 1])
        variants = SeparateMultiallelicVariants.out.map(it -> it[0, 1])
}

workflow CrosscheckContaminantVariants{
    take:
        variants
        synthetics

    main:
        IntersectVCF(variants.combine(synthetics.map(it -> it[1])))

    emit:
        IntersectVCF.out.map(it -> it[0, 1])
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
        // AnnotateVCF(FilterFreebayesVariants.out)

    emit:
        locations = Freebayes.out
        variants = FilterFreebayesVariants.out
}

workflow CallVariantsMutect2{
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


workflow CallVariantsValidate{
    // Allow to determine VAF for given genomic positions in both bed and vcf format
    take:
        reads_aligned
        positions
        reference

    main:
        // Determine noisy SNPs, indels
        FindVariants(reference, reads_aligned, positions)
        noisy_snp = FindVariants.out.map(it -> it[0, 1])
        noisy_indel = FindVariants.out.map(it -> it[0, 2])
        SortNoisySnp(noisy_snp)
        SortNoisyIndel(noisy_indel)
        MergeNoisyVCF(SortNoisySnp.out.combine(SortNoisyIndel.out, by: 0), 'noisy')

        // Filter SNPs, indels
        PerbaseBaseDepthConsensus(reads_aligned.combine(reference), positions, 'consensus.tsv')
        FilterSnp(noisy_snp.combine(PerbaseBaseDepthConsensus.out, by: 0), 'snp')
        FilterIndel(noisy_indel.combine(PerbaseBaseDepthConsensus.out, by: 0), 'indel')
        SortFilteredSnp(FilterSnp.out)
        SortFilteredIndel(FilterIndel.out)
        MergeFilteredVCF(SortFilteredSnp.out.combine(SortFilteredIndel.out, by: 0), 'filtered')

        // Annotate VCF
        AnnotateVCF(MergeFilteredVCF.out)

    emit:
        locations = MergeNoisyVCF.out.map(it -> it[0, 1])
        variants = AnnotateVCF.out.map(it -> it[0, 1])
}
