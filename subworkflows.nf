nextflow.enable.dsl = 2

include {
    // Reference, indexing and input parsing
    Reference_info
    MergeFasta
    BwaIndex
    Minimap2Index as IndexReference
    Minimap2Index as IndexCombined
    FilterShortReads
    SplitReadsOnAdapterSequence
    SplitReadFilesOnNumberOfReads

    // Alignment, annotation
    Minimap2AlignAdaptiveParameterized
    Minimap2Align
    BwaMemSorted
    BwaMemContaminants
    AnnotateBamXTags
    AnnotateBamYTags
    BamAlignmentRateFilter
    PrimaryMappedFilter
    SamtoolsIndexWithID
    MapqAndNMFilter
    BamTagFilter
    FindRegionOfInterest

    // Consensus
    Cycas

    // Variant calling and contaminantion QC
    IntersectVCF
    FindVariants
    SortVCF as SortNoisySnp
    SortVCF as SortNoisyIndel
    SortVCF as SortFilteredSnp
    SortVCF as SortFilteredIndel
    MergeVCF as MergeNoisyVCF
    MergeVCF as MergeFilteredVCF
    FilterValidateVariants as FilterSnp
    FilterValidateVariants as FilterIndel
    AnnotateVCF
    FreebayesContaminants
    SeparateMultiallelicVariants

    // Reporting
    CountFastqInfo as FastqInfoRaw
    CountFastqInfo as FastqInfoConsensus
    PlotFastqsQUalAndLength as PlotRawFastqHist
    PlotFastqsQUalAndLength as PlotFilteredHist
    PlotFastqsQUalAndLength as PlotConFastqHist
    PlotReadStructure
    SamtoolsQuickcheck
    SamtoolsIdxStats
    CountNonBackboneVariants
    PlotMetadataStats
    PerbaseBaseDepth as PerbaseBaseDepthSplit
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
    PlotQScores
    PlotVcf
    SamtoolsMergeBams as SamtoolsMergeBams
    SamtoolsMergeBams as SamtoolsMergeBamsFiltered
    SamtoolsFlagstats
    PasteVariantTable
    PasteContaminantTable
    PlotReport
} from "./processes.nf"

workflow PrepareGenome {
    take:
        reference_genome
        reference_genome_name
        backbones_fasta
    main:
        if (reference_genome_name.endsWith('.txt')) {
            println("txt reference, not implemented. Exiting...")
            genome = "missing"
            exit(1)
        }
        else if (reference_genome_name.endsWith('.gz')) {
            println("gzipped reference, not implemented. Exiting...")
            genome = "missing"
            exit(1)
        }
        else {
            genome = reference_genome
        }
        MergeFasta(genome, backbones_fasta)
        IndexCombined(MergeFasta.out)
        Reference_info(MergeFasta.out)
        
    emit:
        mmi_combi = IndexCombined.out
        fasta_combi = MergeFasta.out
        fasta_ref = reference_genome
}

workflow FilterWithAdapterDetection {
    take:
    read_fq_ch

    main:
    if (params.split_on_adapter == true) {
        fastq = SplitReadsOnAdapterSequence(read_fq_ch)
    }
    else {
        fastq = read_fq_ch
    }

    FilterShortReads(fastq)

    if (params.split_fastq_by_size == true) {
        filtered_reads = SplitReadFilesOnNumberOfReads(FilterShortReads.out)
        // The path to split file is now a list of files,
        // We need to explode the result such that each element is complete
        // With a sample id, file id and path to split file

        filtered_reads = filtered_reads.flatMap { sample_id, file_id, file_paths ->
            file_paths.collect { file_path ->
                def new_file_id = file_path.getBaseName()
                // ideally new files would not end in .part_001.fastq, but instead _part_001.fastq
                // Need new seqkit version to fix this
                return [sample_id, new_file_id, file_path]
            }
        }

    } else {
        filtered_reads = FilterShortReads.out
    }

    emit:
    filtered_reads
}

workflow CycasConsensus{
    take:
        read_fastq
        reference_genome

    main:
        Minimap2AlignAdaptiveParameterized(read_fastq, reference_genome)
        SamtoolsIndexWithID(Minimap2AlignAdaptiveParameterized.out)
        PrimaryMappedFilter(SamtoolsIndexWithID.out) 
        MapqAndNMFilter(PrimaryMappedFilter.out)
        Cycas(MapqAndNMFilter.out)

    emit:
        fastq = Cycas.out.map(it -> tuple(it[0], it[1], it[2]))
        json = Cycas.out.map(it -> tuple(it[0], it[1], it[3]))
        split_bam = Minimap2AlignAdaptiveParameterized.out
        split_bam_filtered = MapqAndNMFilter.out.map(it -> tuple(it[0], it[1], it[2]))
}

workflow AlignByID {
    take:
        reads
        reference_genome
        jsons

    main:
        id = reads.map(it -> it[0])
        Minimap2Align(reads, reference_genome)
        metadata_pairs = Minimap2Align.out.join(jsons).map(it -> tuple(it[0], it[1], it[2], it[4]))
	    AnnotateBamYTags(metadata_pairs)
        bams = AnnotateBamYTags.out.groupTuple(by: 0).map(it -> tuple(it[0], it[1], it[2]))
        bam = SamtoolsMergeBams(bams)

    emit:
        bam
}

workflow AlignContaminants {
    take:
        synthetic_reads
        reference_genome
        reference_genome_indexes

    main:
        bwa_index_file_count = 5
        // We do a smaller than since there might be a .fai file as well!
        if (reference_genome_indexes.size < bwa_index_file_count){
            println "==================================="
            println "Warning! BWA index files are missing for the reference genome, This will slowdown execution in a major way."
            println "==================================="
            println ""
            println ""
            reference_genome_indexes = BwaIndex(reference_genome)
        }
        BwaMemSorted(synthetic_reads, reference_genome, reference_genome_indexes.collect() )

    emit:
        bam = BwaMemSorted.out
}

workflow AnnotateBam {
    take:
        reads
        sequencing_summary

    main:
        AnnotateBamXTags(reads, sequencing_summary)
    emit:
        AnnotateBamXTags.out
}

workflow FilterBam {
    take:
        annotated_bam
        minimun_repeat_count

    main:
        BamTagFilter(annotated_bam, 'YM', minimun_repeat_count)
        if (params.region_file == 'auto') {
            filtered_bam = BamTagFilter.out
        } else {
            BamAlignmentRateFilter(BamTagFilter.out)
            filtered_bam = BamAlignmentRateFilter.out
        }
    emit:
        filtered_bam
}

// workflow ProcessTargetRegions{
//     take:
//         reads_aligned
    
//     main:
//         positions = FindRegionOfInterest(reads_aligned)

//     emit:
//         positions

// }

workflow ProcessTargetRegions{
    take:
        variant_file_name
        reads_aligned
    
    main:
        variant_file = Channel.fromPath(params.region_file, checkIfExists: false)

        if (variant_file_name.endsWith('.bed')) {
            positions = variant_file
        } else {
            // 'auto' or otherwise
            positions = FindRegionOfInterest(reads_aligned)
        }

    emit:
        positions
}

workflow CallVariantsValidate{
    // Allow to determine VAF for given genomic positions in both bed and vcf format
    take:
        reads_aligned
        positions
        reference

    main:
        // Determine noisy SNPs, indels
        FindVariants(reference, reads_aligned.combine(positions, by: [0, 1]))
        noisy_snp = FindVariants.out.map(it -> it[0, 1, 2])
        noisy_indel = FindVariants.out.map(it -> it[0, 1, 3])
        SortNoisySnp(noisy_snp)
        SortNoisyIndel(noisy_indel)
        MergeNoisyVCF(SortNoisySnp.out.combine(SortNoisyIndel.out, by: [0, 1]), 'noisy') // TODO: verify that sample_id == file_id

        // Filter SNPs, indels
        PerbaseBaseDepthConsensus(reads_aligned.combine(reference), positions, 'consensus.tsv')
        FilterSnp(noisy_snp.combine(PerbaseBaseDepthConsensus.out, by: [0, 1]), 'snp')
        FilterIndel(noisy_indel.combine(PerbaseBaseDepthConsensus.out, by: [0, 1]), 'indel')
        SortFilteredSnp(FilterSnp.out)
        SortFilteredIndel(FilterIndel.out)
        MergeFilteredVCF(SortFilteredSnp.out.combine(SortFilteredIndel.out, by: [0, 1]), 'filtered')

        // Annotate VCF
        AnnotateVCF(MergeFilteredVCF.out)

    emit:
        locations = MergeNoisyVCF.out
        variants = AnnotateVCF.out
}

workflow CallContaminantMutants{
    take:
        reads_aligned
        positions
        reference

    main:
        indexed_ref = IndexReference(reference)
        FreebayesContaminants(reads_aligned.combine(indexed_ref), positions)
        SeparateMultiallelicVariants(FreebayesContaminants.out)

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

workflow Report {
    take:
        reference_fasta
        fastq_raw
        fastq_filtered
        split_bam
        split_bam_filtered
        fastq_consensus
        read_info
        consensus_bam
        roi
        noisy_vcf
        variants_vcf
        contaminants_vcf

    main:
        // Raw FASTQ reads
        fastq_raw_extension = fastq_raw.map(it -> it[2]).first().getExtension()
        fastq_raw_grouped = fastq_raw.map(it -> it[0, 2]).groupTuple(by: 0)
        FastqInfoRaw(fastq_raw_grouped, 'raw')
        PlotRawFastqHist(fastq_raw_grouped, fastq_raw_extension, "raw", '"Raw fastq info"')
        
        // Fitletered FASTQ reads
        fastq_filtered_extension = fastq_filtered.map(it -> it[2]).first().getExtension()
        fastq_filtered_grouped = fastq_filtered.map(it -> it[0, 2]).groupTuple(by: 0)
        PlotFilteredHist(fastq_filtered_grouped, fastq_filtered_extension, "filtered", '"Filtered fastq info"')

        // Consensus FASTQ reads
        fastq_consensus_extension = fastq_consensus.map(it -> it[2]).first().getExtension()
        fastq_consensus_grouped = fastq_consensus.map(it -> it[0, 2]).groupTuple(by: 0)
        FastqInfoConsensus(fastq_consensus_grouped, 'consensus')
        PlotConFastqHist(fastq_consensus_grouped, fastq_consensus_extension, "consensus", '"Consensus fastq info"')

        // Split BAM
        merged_split_bam = SamtoolsMergeBams(split_bam)
        PlotReadStructure(merged_split_bam) // 90
        
        // Consensus BAM
        SamtoolsQuickcheck(consensus_bam)
        SamtoolsIdxStats(consensus_bam)
        CountNonBackboneVariants(variants_vcf)

        // Metadata
        meta_data = read_info.map(it -> it[0, 2]).groupTuple(by: 0)
        PlotMetadataStats(meta_data) // 92 // TODO: output per sample

        // roi = FindRegionOfInterest(consensus_bam)
      
        // QScores
        PerbaseBaseDepthSplit(merged_split_bam.combine(reference_fasta), roi, 'split.tsv')
        PerbaseBaseDepthConsensus(consensus_bam.combine(reference_fasta), roi, 'consensus.tsv')
        qscores = PerbaseBaseDepthSplit.out.combine(PerbaseBaseDepthConsensus.out, by: [0, 1]).map(it -> it[0, 2, 3])
        PlotQScores(qscores) // 91
        
        // VCF
        PlotVcf(noisy_vcf)

        // Filtered Split BAM
        merged_split_bam_filtered = SamtoolsMergeBamsFiltered(split_bam_filtered)
        SamtoolsFlagstats(merged_split_bam_filtered)

        // Variants
        PasteVariantTable(variants_vcf)

        // Contaminants
        // PasteContaminantTable(contaminants_vcf)

        // Plot Report
        PlotReport(
            PlotRawFastqHist.out.combine(
            PlotFilteredHist.out, by: 0).combine(
            PlotConFastqHist.out, by: 0).combine(
            PlotReadStructure.out, by: 0).combine(
            PlotQScores.out, by: 0).combine(
            PlotVcf.out, by: 0).combine(
            PlotMetadataStats.out, by: 0).combine(
            PasteVariantTable.out, by: 0).combine(
            // PasteContaminantTable.out, by: 0).combine(
            SamtoolsFlagstats.out, by: 0).combine(
            CountNonBackboneVariants.out, by: 0).combine(
            SamtoolsIdxStats.out, by: 0
            )
        )
}