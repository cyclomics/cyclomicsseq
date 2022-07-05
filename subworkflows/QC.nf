
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
} from "./modules/bin"

include {
    CountFastqInfo as FastqInfoRaw
    CountFastqInfo as FastqInfoConsensus
} from "./modules/seqkit"

include {
    SamtoolsQuickcheck
    SamtoolsFlagstats
} from "./modules/samtools"

include {
    FindRegionOfInterest
} from "./modules/mosdepth"

include {
    SamtoolsMergeBams
} from "./modules/samtools"

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
        read_info
        fastq_raw
        fastq_consensus
        consensus_bam
        quick_results
        split_bam

    main:
        // dont add the ID to the process
        // CollectClassificationTypes(read_info.map(it -> it[1]).collect())
        FastqInfoRaw(fastq_raw.collect())
        PlotRawFastqHist(fastq_raw.collect(), "*.fastq.gz")
        FastqInfoConsensus(fastq_consensus.map(it -> it[1]).collect())
        PlotConFastqHist(fastq_consensus.map(it -> it[1]).collect(), "*.fastq")
        
        SamtoolsMergeBams('splibams_merged', split_bam.collect())
        PlotReadStructure(SamtoolsMergeBams.out)

        SamtoolsQuickcheck(consensus_bam)
        SamtoolsFlagstats(consensus_bam)
        FindRegionOfInterest(consensus_bam)

        if (quick_results == true) {
            SamtoolsQuickcheck.out.view()
            SamtoolsFlagstats.out.view()
        }
}
