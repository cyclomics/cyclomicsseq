
nextflow.enable.dsl=2


include {
    BwaMemSorted
} from "../modules/bwa.nf"

include {
    LastCreateDB
    LastTrainModelFastq
    LastTrainModelFasta
    LastAlign
    LastAlignTrained
    LastSplit
    SamtoolsFixSam
} from "../modules/last.nf"

include{
    ConcatenateFasta
    ConcatenateFastq
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

workflow LastalAlignFasta{
    take:
        read_fq_ch
        reference_genome
    main:
        ConcatenateFasta(read_fq_ch.collect())
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

workflow LastalAlignFastq{
    take:
        read_fq_ch
        reference_genome
    main:
        ConcatenateFastq(read_fq_ch.collect())
        LastCreateDB(reference_genome)
        // CreateDB makes many db.* files, need all of them downstream
        last_db_collection = LastCreateDB.out.collect()
        // pass the read_fq into lastal
        LastAlign(ConcatenateFastq.out, 
            last_db_collection,
            )
    emit:
        LastAlign.out
}

workflow LastalAlignTrainedFasta{
    take:
        read_fq_ch
        reference_genome
    main:
        LastCreateDB(reference_genome)
        // CreateDB makes many db.* files, need all of them downstream
        last_db_collection = LastCreateDB.out.collect()
        
        reads = ConcatenateFasta(read_fq_ch.collect())
        model =  LastTrainModelFasta(read_fq_ch, last_db_collection)

        // pass the read_fq into lastal
        alignment = LastAlignTrained(reads, 
            last_db_collection,
            model
            )
        SamtoolsFixSam(alignment, reference_genome)
        SamtoolsIndex(SamtoolsFixSam.out)

    emit:
        SamtoolsIndex.out
}

workflow LastalAlignTrainedFastq{
    take:
        read_fq_ch
        reference_genome
    main:
        LastCreateDB(reference_genome)
        // CreateDB makes many db.* files, need all of them downstream
        last_db_collection = LastCreateDB.out.collect()
        reads = ConcatenateFastq(read_fq_ch.collect())
        model =  LastTrainModelFastq(read_fq_ch, last_db_collection)
        alignment = LastAlignTrained(reads, 
            last_db_collection,
            model
            )
        SamtoolsFixSam(alignment, reference_genome)
        SamtoolsIndex(SamtoolsFixSam.out)

    emit:
        SamtoolsIndex.out
}