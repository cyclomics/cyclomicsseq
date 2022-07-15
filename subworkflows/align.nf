
nextflow.enable.dsl=2

include {
    AnnotateBamXTags
    AnnotateBamYTags
} from "./modules/bin.nf"

include {
    BwaMemSorted
} from "./modules/bwa.nf"

include {
    MinimapAlign
    Minimap2AlignAdaptive
    Minimap2Index as IndexReference
    Minimap2Index as IndexCombined
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
    SamtoolsDepthToTSV
    SamtoolsMergeBams
} from "./modules/samtools"

include {
    MergeFasta
} from "./modules/seqkit"


workflow PrepareGenome {
    take:
        reference_genome
        reference_genome_name
        backbones_fasta
    main:
        if (reference_genome_name.endsWith('.txt')) {
            println("txt reference, not implemented. Exiting...")
            genome = "missing"
            exit(1)
        }
        else if (reference_genome_name.endsWith('.gz')) {
            println("gzipped reference, not implemented. Exiting...")
            genome = "missing"
            exit(1)
        }
        else {
            genome = reference_genome
        }
        MergeFasta(genome, backbones_fasta)
        IndexCombined(MergeFasta.out)
        IndexReference(genome)
    emit:
        mmi_combi = IndexCombined.out
        mmi_ref = IndexReference.out
        fasta_combi = MergeFasta.out
        fasta_ref = reference_genome
}

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
        jsons

    main:

        Minimap2AlignAdaptive(reads.map(it -> it[1]), reference_genome)
        id = reads.first()map( it -> it[0])
        id = id.map(it -> it.split('_')[0])

        AnnotateBamYTags(Minimap2AlignAdaptive.out.join(jsons))
        bams = Minimap2AlignAdaptive.out.map(it -> it[1]).collect()
        
        SamtoolsMergeBams(id, bams)

    emit:
        bam = SamtoolsMergeBams.out
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

workflow Annotate{
    take:
        reads
        sequencing_summary
    main:
        AnnotateBamXTags(reads, sequencing_summary)
    emit:
        AnnotateBamXTags.out
}