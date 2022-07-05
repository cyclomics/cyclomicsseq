
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
        reference_fasta
        fastq_raw
        split_bam
        fastq_consensus
        read_info
        consensus_bam
        quick_results

    main:
        // dont add the ID to the process
        // CollectClassificationTypes(read_info.map(it -> it[1]).collect())
        first_fq = fastq_raw.first()
        id = first_fq.simpleName
        extension = first_fq.getExtension()

        println(id)
        println(extension)
        FastqInfoRaw(fastq_raw.collect())
        PlotRawFastqHist(fastq_raw.collect(), extension, id)
        
        id = fastq_consensus.first().map(it -> it[0])
        extension = fastq_consensus.map(it -> it[1]).first().getExtension()

        FastqInfoConsensus(fastq_consensus.map(it -> it[1]).collect())
        PlotConFastqHist(fastq_consensus.map(it -> it[1]).collect(),extension, id)

        SamtoolsMergeBams('splibams_merged', split_bam.collect())
        PlotReadStructure(SamtoolsMergeBams.out)
        
        SamtoolsQuickcheck(consensus_bam)
        SamtoolsFlagstats(consensus_bam)
        FindRegionOfInterest(consensus_bam)
        // if (params.variant_calling == "validate") {
        //     PlotVcf(vcf)
        // }
        if (quick_results == true) {
            SamtoolsQuickcheck.out.view()
            SamtoolsFlagstats.out.view()
        }
}
