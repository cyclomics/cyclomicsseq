
nextflow.enable.dsl=2


include {DummyProcess as Consensus;
} from "../modules/dummy.nf"


workflow  ConsensusBasic{
    take:
        read_fq_ch
        backbone_fq_ch
    emit:
        Consensus.output
    main:
        Consensus(read_fq_ch)
        
}
