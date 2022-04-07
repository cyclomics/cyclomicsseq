
nextflow.enable.dsl=2


include {
    BwaMemSorted
} from "./modules/bwa.nf"

include {
    MinimapAlign
} from "./modules/minimap.nf"

include {
    LastCreateDB
    LastTrainModelFastq
    LastTrainModelFasta
    LastAlign
    LastAlignTrained
    LastSplit
    Maf2sam
    SamtoolsFixSam
} from "./modules/last.nf"

include{
    ConcatenateFasta
    ConcatenateFastq
} from "./modules/utils.nf"

include {
    SamtoolsIndex
    SamToBam
    SamtoolsMergeTuple
    SamtoolsDepth
    SamtoolsDepthToJson
    SamtoolsMergeBams
} from "./modules/samtools"


workflow AlignBWA{
    take:
        read_fq_ch
        reference_genome
    main:
        BwaMemSorted(read_fq_ch, reference_genome)
    emit:
        BwaMemSorted.out
}

workflow Minimap2Align{
    // Call minimap2 on all reads files (tuple(x,bam)) convert to bam and merge using samtools
    take:
        reads
        reference_genome
    main:
        MinimapAlign(reads.combine(reference_genome))
        SamToBam(MinimapAlign.out)
        id = reads.first()map( it -> it[0])
        id = id.map(it -> it.split('_')[0])
        bams = SamToBam.out.map(it -> it[1]).collect()
        SamtoolsMergeBams(id, bams)
        SamtoolsDepth(SamtoolsMergeBams.out)
        SamtoolsDepthToJson(SamtoolsDepth.out)

    emit:
        bam = SamtoolsMergeBams.out
        depth = SamtoolsDepthToJson.out
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
        alignment = LastAlign(ConcatenateFasta.out, 
            last_db_collection,
            )
        samfile = Maf2sam(alignment)
        SamtoolsFixSam(samfile, reference_genome)
        SamtoolsIndex(SamtoolsFixSam.out)
        
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
        alignment = LastAlign(ConcatenateFastq.out, 
            last_db_collection,
        )
            samfile = Maf2sam(alignment)
            SamtoolsFixSam(samfile, reference_genome)
            SamtoolsIndex(SamtoolsFixSam.out)

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
        model =  LastTrainModelFasta(reads, last_db_collection)

        // pass the read_fq into lastal
        alignment = LastAlignTrained(reads, 
            last_db_collection,
            model
            )
        samfile = Maf2sam(alignment)
        SamtoolsFixSam(samfile, reference_genome)
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
        model =  LastTrainModelFastq(reads, last_db_collection)
        alignment = LastAlignTrained(reads, 
            last_db_collection,
            model
            )
        samfile = Maf2sam(alignment)
        SamtoolsFixSam(samfile, reference_genome)
        SamtoolsIndex(SamtoolsFixSam.out)

    emit:
        SamtoolsIndex.out
}