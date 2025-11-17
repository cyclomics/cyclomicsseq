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

    // Reporting
    CountFastqInfo as FastqInfoRaw
    CountFastqInfo as FastqInfoConsensus
    PlotFastqsQUalAndLength as PlotRawFastqHist
    PlotFastqsQUalAndLength as PlotFilteredHist
    PlotFastqsQUalAndLength as PlotConFastqHist
    PlotReadStructure
    SamtoolsQuickcheck
    SamtoolsIdxStats
    PlotMetadataStats
    PerbaseBaseDepth as PerbaseBaseDepthSplit
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
    PlotQScores
    SamtoolsMergeBams as SamtoolsMergeBams
    SamtoolsMergeBams as SamtoolsMergeBamsFiltered
    SamtoolsFlagstats
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

    if (params.split_fastq_by_size == true) {
        split_reads = SplitReadFilesOnNumberOfReads(fastq)

        split_reads = split_reads.flatMap { sample_id, file_id, file_paths ->
            if (file_paths instanceof List) {
                file_paths.collect { file_path ->
                    def new_file_id = file_path.getBaseName()
                    return [sample_id, new_file_id, file_path]
                }
            } else {
                def new_file_id = file_paths.getBaseName()
                return [[sample_id, new_file_id, file_paths]]
            }
        }
    } else {
        split_reads = fastq
    }

    filtered_reads = FilterShortReads(split_reads)

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
        Minimap2Align(reads, reference_genome)
        metadata_pairs = Minimap2Align.out.combine(jsons, by:[0, 1])
	    AnnotateBamYTags(metadata_pairs)
        bams = AnnotateBamYTags.out.groupTuple(by: 0).map(it -> tuple(it[0], it[1], it[2]))
        SamtoolsMergeBams(bams)
        bam = SamtoolsMergeBams.out

    emit:
        bam
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
        minimum_repeat_count

    main:
        BamTagFilter(annotated_bam, 'YM', minimum_repeat_count)
        if (params.region_file != 'auto' && params.filter_by_alignment_rate == true) {
            BamAlignmentRateFilter(BamTagFilter.out)
            filtered_bam = BamAlignmentRateFilter.out
        } else {
            filtered_bam = BamTagFilter.out
        }
    emit:
        filtered_bam
}

workflow ProcessTargetRegions{
    take:
        region_file
        reads_aligned
    
    main:
        regions = Channel.fromPath(region_file, checkIfExists: false)
        positions = FindRegionOfInterest(reads_aligned.combine(regions))

    emit:
        positions
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

    main:
        // Raw FASTQ reads
        fastq_raw_extension = fastq_raw.map(it -> it[2]).first().getExtension()
        fastq_raw_grouped = fastq_raw.map(it -> it[0, 2]).groupTuple(by: 0)
        FastqInfoRaw(fastq_raw_grouped, 'raw')
        PlotRawFastqHist(fastq_raw_grouped, fastq_raw_extension, "raw", '"Raw read info"')
        
        // Fitltered FASTQ reads
        fastq_filtered_extension = fastq_filtered.map(it -> it[2]).first().getExtension()
        fastq_filtered_grouped = fastq_filtered.map(it -> it[0, 2]).groupTuple(by: 0)
        PlotFilteredHist(fastq_filtered_grouped, fastq_filtered_extension, "filtered", '"Filtered read info"')

        // Consensus FASTQ reads
        fastq_consensus_extension = fastq_consensus.map(it -> it[2]).first().getExtension()
        fastq_consensus_grouped = fastq_consensus.map(it -> it[0, 2]).groupTuple(by: 0)
        FastqInfoConsensus(fastq_consensus_grouped, 'consensus')
        PlotConFastqHist(fastq_consensus_grouped, fastq_consensus_extension, "consensus", '"Consensus read info"')

        // Split BAM
        split_bam_grouped = split_bam.groupTuple(by: 0)
        merged_split_bam = SamtoolsMergeBams(split_bam_grouped)
        PlotReadStructure(merged_split_bam) // 90
        
        // Consensus BAM
        SamtoolsQuickcheck(consensus_bam)
        SamtoolsIdxStats(consensus_bam)

        // Metadata
        meta_data = read_info.map(it -> it[0, 2]).groupTuple(by: 0)
        PlotMetadataStats(meta_data) // 92
      
        // QScores
        PerbaseBaseDepthSplit(merged_split_bam.combine(roi, by: [0, 1]), reference_fasta, 'split.tsv')
        PerbaseBaseDepthConsensus(consensus_bam.combine(roi, by: [0, 1]), reference_fasta, 'consensus.tsv')
        qscores = PerbaseBaseDepthSplit.out.combine(PerbaseBaseDepthConsensus.out, by: [0, 1]).map(it -> it[0, 2, 3])
        PlotQScores(qscores) // 91

        // Filtered Split BAM
        split_bam_filtered_grouped = split_bam_filtered.groupTuple(by: 0)
        merged_split_bam_filtered = SamtoolsMergeBamsFiltered(split_bam_filtered_grouped)
        SamtoolsFlagstats(merged_split_bam_filtered)

        // Plot Report
        PlotReport(
            PlotRawFastqHist.out
            .combine(PlotFilteredHist.out, by: 0)
            .combine(PlotConFastqHist.out, by: 0)
            .combine(PlotReadStructure.out, by: 0)
            .combine(PlotQScores.out, by: 0)
            .combine(PlotMetadataStats.out, by: 0)
            .combine(SamtoolsFlagstats.out, by: 0)
            .combine(SamtoolsIdxStats.out, by: 0)
            .map { values ->
                def key = values[0]
                def files = values[1..-1]
                return [key, files]
            }
        )
}