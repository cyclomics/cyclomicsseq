
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
} from "./modules/utils"

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
    main:
        CollectClassificationTypes(read_info.collect())
        // CountReadsRaw()
        // CountReadsConsensus()
}