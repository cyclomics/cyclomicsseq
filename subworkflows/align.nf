
nextflow.enable.dsl=2


include {
    BwaMemSorted
} from "../modules/bwa.nf"

include {
    LastCreateDB
    LastTrainModel
    LastAlign
    LastSplit
    SamtoolsFixSam
} from "../modules/last.nf"

include{
    ConcatenateFasta
} from "../modules/utils.nf"

include {
    SamtoolsIndex
} from "../modules/samtools"


workflow AlignBWA{
    take:
        read_fq_ch
        reference_genome
    main:
        BwaMemSorted(read_fq_ch, reference_genome)
    emit:
        BwaMemSorted.out
}

workflow LastalAlign{
    take:
        read_fq_ch
        reference_genome
    main:
        ConcatenateFasta(read_fq_ch)
        LastCreateDB(reference_genome)
        // CreateDB makes many db.* files, need all of them downstream
        last_db_collection = LastCreateDB.out.collect()
        // pass the read_fq into lastal
        LastAlign(ConcatenateFasta.out, 
            last_db_collection,
            )
    emit:
        LastAlign.out
}

workflow LastalAlignTrained{
    take:
        read_fq_ch
        reference_genome
    main:
        ConcatenateFasta(read_fq_ch)
        LastCreateDB(reference_genome)
        // CreateDB makes many db.* files, need all of them downstream
        last_db_collection = LastCreateDB.out.collect()
        LastTrainModel(read_fq_ch, last_db_collection)

        // pass the read_fq into lastal
        LastAlign(ConcatenateFasta.out, 
            last_db_collection,
            )
        Maf2sam(LastAlign.out)
        SamtoolsFixSam(Maf2sam.out, reference_genome)
        SamtoolsIndex(SamtoolsFixSam.out)

    emit:
        SamtoolsIndex.out
}