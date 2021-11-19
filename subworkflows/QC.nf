
include {PYCOQC;
} from "../modules/pycoqc.nf"


workflow  QC_pycoqc{
    take:
        read_directory
    emit:
        PYCOQC.output
        PYCOQC.output
    main:
        PYCOQC(read_directory)
}
