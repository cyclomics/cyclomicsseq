
nextflow.enable.dsl=2

include {
    PYCOQC
} from "./modules/pycoqc.nf"

include {
    MinionQc
    MinionQcToJson
} from "./modules/minion_qc.nf"

include {
    CollectClassificationTypes
    PlotReadStructure
    PlotFastqsQUalAndLength as PlotRawFastqHist
    PlotFastqsQUalAndLength as PlotConFastqHist
    PlotVcf
    PlotQScores
    PlotMetadataStats
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
    SamtoolsMergeBams
} from "./modules/samtools"

include {
    PerbaseBaseDepth as PerbaseBaseDepthSplit
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
} from "./modules/perbase"


workflow  QC_pycoqc{
    take:
        read_directory
    main:
        PYCOQC(read_directory)
    emit:
        PYCOQC.output
}

workflow QC_MinionQc {
    take:
        read_directory
    main:
        MinionQc(read_directory)
        MinionQcToJson(MinionQc.out.summary)
    emit:
        MinionQcToJson.output
}


workflow PostQC {
    take:
        reference_fasta
        fastq_raw
        split_bam
        fastq_consensus
        read_info
        consensus_bam
        quick_results
        vcf

    main:
        // dont add the ID to the process
        // CollectClassificationTypes(read_info.map(it -> it[1]).collect())
        first_fq = fastq_raw.first()
        id = first_fq.simpleName
        extension = first_fq.getExtension()

        FastqInfoRaw(fastq_raw.collect())
        PlotRawFastqHist(fastq_raw.collect(), extension, id)
        
        first_fq = fastq_consensus.first()
        id = first_fq.map(it -> it[0])
        extension = first_fq.map(it -> it[1]).getExtension()

        FastqInfoConsensus(fastq_consensus.map(it -> it[1]).collect())
        PlotConFastqHist(fastq_consensus.map(it -> it[1]).collect(),extension, id)

        merged_split_bam = SamtoolsMergeBams('splibams_merged', split_bam.collect())
        PlotReadStructure(merged_split_bam)
        PlotMetadataStats(id, read_info.map(it -> it[1]).collect())

        roi = FindRegionOfInterest(consensus_bam)
        // pileups_split_bam = MPileupSplit(merged_split_bam.combine(reference_fasta), roi)
        // pileups_consensus_bam = MPileupConsensus( consensus_bam.combine(reference_fasta), roi)
        
        PerbaseBaseDepthSplit(merged_split_bam.combine(reference_fasta), roi, 'split.tsv')
        PerbaseBaseDepthConsensus(consensus_bam.combine(reference_fasta), roi, 'consensus.tsv')

        PlotQScores(PerbaseBaseDepthSplit.out, PerbaseBaseDepthConsensus.out)

        if (params.variant_calling == "validate") {
            PlotVcf(vcf)
        }
        
        SamtoolsQuickcheck(consensus_bam)
        SamtoolsFlagstats(consensus_bam)
        
        if (quick_results == true) {
            SamtoolsQuickcheck.out.view()
            SamtoolsFlagstats.out.view()
        }
}


// workflow CreateAccuracyPlot {
//     take:
//         merged_consensus_bam
//         merged_split_bam
//         reference_fasta

//     main:


// }