
nextflow.enable.dsl=2


include {DummyProcess as Repeat_ID;
} from "../modules/dummy.nf"

include {
    Tidehunter
} from "../modules/tidehunter.nf"

workflow  TideHunterBasic{
    take:
        read_fq_ch
    emit:
        Tidehunter
    main:
        Tidehunter(read_fq_ch)
}
