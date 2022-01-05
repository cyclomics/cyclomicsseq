
include {PYCOQC;
} from "../modules/pycoqc.nf"

include {MinionQc;
} from "../modules/minion_qc.nf"


workflow  QC_pycoqc{
    take:
        read_directory
    emit:
        PYCOQC.output
    main:
        PYCOQC(read_directory)
}

workflow QC_MinionQc {
    take:
        read_directory
    emit:
        MinionQc.output
    main:
        MinionQc(read_directory)
}