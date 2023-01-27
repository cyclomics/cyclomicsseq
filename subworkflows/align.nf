
nextflow.enable.dsl=2

include {
    AnnotateBamXTags
    AnnotateBamYTags
} from "./modules/bin.nf"

include {
    BwaIndex
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
    BamTagFilter
} from "./modules/samtools"

include {
    MergeFasta
} from "./modules/seqkit"



process Reference_info {
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta

    output:
        path "reference_info.txt"

    script:
        """
        echo "--------- md5 hash info ---------" >> reference_info.txt
        md5sum $fasta >> reference_info.txt
        echo "--------- assembly info ---------" >> reference_info.txt
        seqkit fx2tab --length --name --header-line --seq-hash $fasta >> reference_info.txt
        """
}

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
        // IndexReference(genome)
        Reference_info(MergeFasta.out)
        
    emit:
        mmi_combi = IndexCombined.out
        // mmi_ref = IndexReference.out
        fasta_combi = MergeFasta.out
        fasta_ref = reference_genome
}

workflow BWAAlign{
    take:
        reads
        reference_genome
        reference_genome_indexes
        jsons

    main:
        bwa_index_file_count = 5
        // We do a smaller than since there might be a .fai file as well!
        if (reference_genome_indexes.size < bwa_index_file_count){
            println "Warning! BWA index files are missing for the reference genome, This will slowdown execution in a major way."
            reference_genome_indexes = BwaIndex(reference_genome)
        }
        
        id = reads.first().map( it -> it[0])
        id = id.map(it -> it.split('_')[0])
        reads_fastq = reads.map(it -> it[1])
        BwaMemSorted(reads_fastq, reference_genome, reference_genome_indexes.collect() )

        if (params.consensus_calling == "cycas") {
            // For now we only do Y tag addition for cycas
            metadata_pairs = BwaMemSorted.out.join(jsons)
	        AnnotateBamYTags(metadata_pairs)
            bams = AnnotateBamYTags.out.map(it -> it[1]).collect()
        }
        else {
            bams= Minimap2AlignAdaptive.out
        }

        SamtoolsMergeBams(id, bams.collect())

    emit:
        bam = SamtoolsMergeBams.out
}

workflow Minimap2Align{
    // Call minimap2 on all reads files (tuple(x,bam)) convert to bam and merge using samtools.
    // Adds metadata tags (eg YM) from metadata
    // Will create a single bam given filenames like FAS12345_blabla : FAS12345.bam
    take:
        reads
        reference_genome
        jsons
        consensus_calling

    main:
        Minimap2AlignAdaptive(reads.map(it -> it[1]), reference_genome)
        id = reads.first().map( it -> it[0])
        id = id.map(it -> it.split('_')[0])

        if (consensus_calling == "cycas") {
            // For now we only do Y tag addition for cycas
            metadata_pairs = Minimap2AlignAdaptive.out.join(jsons)
            AnnotateBamYTags(metadata_pairs)
            bams = AnnotateBamYTags.out.map(it -> it[1]).collect()
        }
        else {
            bams= Minimap2AlignAdaptive.out
        }
        
        SamtoolsMergeBams(id, bams)

    emit:
        bam = SamtoolsMergeBams.out
}

workflow AnnotateFilter{
    take:
        reads
        sequencing_summary
        minimun_repeat_count

    main:
        AnnotateBamXTags(reads, sequencing_summary)
        BamTagFilter(AnnotateBamXTags.out, 'YM', minimun_repeat_count)
    emit:
        BamTagFilter.out
}
