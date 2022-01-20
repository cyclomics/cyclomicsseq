
nextflow.enable.dsl=2


include {DummyProcess as Filtering;
} from "./modules/dummy.nf"

workflow  FilteringBasic{
    take:
        read_fq_ch
    emit:
        Filtering.output
    main:
        Filtering(read_fq_ch)
}
