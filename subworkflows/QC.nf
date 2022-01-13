
nextflow.enable.dsl=2

include {
    PYCOQC
} from "../modules/pycoqc.nf"

include {
    MinionQc
    MinionQcToJson
} from "../modules/minion_qc.nf"


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
