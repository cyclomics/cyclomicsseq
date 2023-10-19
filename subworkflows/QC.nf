
nextflow.enable.dsl=2

include {
    CollectClassificationTypes
    PlotReadStructure
    PlotFastqsQUalAndLength as PlotRawFastqHist
    PlotFastqsQUalAndLength as PlotFilteredHist
    PlotFastqsQUalAndLength as PlotConFastqHist
    PlotVcf
    PasteVariantTable
    PlotQScores
    PlotMetadataStats
    PlotReportDetailed
    PlotReportStd
} from "./modules/bin"

include {
    CountFastqInfo as FastqInfoRaw
    CountFastqInfo as FastqInfoConsensus
} from "./modules/seqkit"

include {
    SamtoolsQuickcheck
    SamtoolsFlagstats
    FindRegionOfInterest
    MPileup as MPileupSplit
    MPileup as MPileupConsensus
} from "./modules/samtools"

include {
    SamtoolsMergeBams as SamtoolsMergeBams
    SamtoolsMergeBams as SamtoolsMergeBamsFiltered
    SamtoolsIdxStats
} from "./modules/samtools"

include {
    PerbaseBaseDepth as PerbaseBaseDepthSplit
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
} from "./modules/perbase"

include{
    CountNonBackboneVariants
} from "./modules/utils"


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

    main:
        // dont add the ID to the process
        first_fq = fastq_raw.first()
        id = first_fq.simpleName
        extension = first_fq.getExtension()
        FastqInfoRaw(fastq_raw.collect(),'raw')
        PlotRawFastqHist(fastq_raw.collect(), extension, id + "raw", '"Raw fastq info"')
        
        first_fq = fastq_filtered.first()
        id = first_fq.simpleName
        extension = first_fq.getExtension()
        PlotFilteredHist(fastq_filtered.collect(),extension, id + "filtered", '"Filtered fastq info"')

        first_fq = fastq_consensus.first()
        id = first_fq.map(it -> it[0])
        extension = first_fq.map(it -> it[1]).getExtension()

        FastqInfoConsensus(fastq_consensus.map(it -> it[1]).collect(), 'consensus')
        PlotConFastqHist(fastq_consensus.map(it -> it[1]).collect(),extension, id + "consensus", '"Consensus fastq info"')

        merged_split_bam = SamtoolsMergeBams('splibams_merged', split_bam.collect())
        PlotReadStructure(merged_split_bam) // 90
        
        SamtoolsQuickcheck(consensus_bam)
        SamtoolsIdxStats(consensus_bam)
        CountNonBackboneVariants(variants_vcf)

        meta_data = read_info.map(it -> it[1]).collect()
        PlotMetadataStats(meta_data) // 92

        // roi = FindRegionOfInterest(consensus_bam)
      
        PerbaseBaseDepthSplit(merged_split_bam.combine(reference_fasta), roi, 'split.tsv')
        PerbaseBaseDepthConsensus(consensus_bam.combine(reference_fasta), roi, 'consensus.tsv')
        PlotQScores(PerbaseBaseDepthSplit.out, PerbaseBaseDepthConsensus.out) // 91
        
        PlotVcf(noisy_vcf)

        merged_split_bam_filtered = SamtoolsMergeBamsFiltered('splibams_filtered_merged',split_bam_filtered.collect())
        SamtoolsFlagstats(merged_split_bam_filtered)

        PasteVariantTable(variants_vcf)
        PlotReportDetailed(
            PlotRawFastqHist.out.combine(
            PlotFilteredHist.out).combine(
            PlotConFastqHist.out).combine(
            PlotReadStructure.out).combine(
            PlotQScores.out).combine(
            PlotVcf.out).combine(
            PlotMetadataStats.out).combine(
            PasteVariantTable.out).combine(
            SamtoolsFlagstats.out).combine(
            CountNonBackboneVariants.out).combine(
            SamtoolsIdxStats.out
            )
        )
}
