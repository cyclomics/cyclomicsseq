#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/CycloSeq
========================================================================================
    Github : https://github.com/cyclomics/cycloseq
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER VALUES
========================================================================================
*/
// ### PARAMETERS
params.input_reads           = "/home/dami/data/raw_data/JAG6252/fastq_pass/FAO43712_pass_75edb72a_?.fastq"

params.backbone_locus        = "BB22"
params.backbone_fasta        = "/home/dami/data/backbones/backbones/BB22.fa"
params.backbone_prime_length = 40

params.insert_target         = "TP53:1-27760"
params.insert_target_genome  = "/home/dami/data/ref_genomes/hg_19/tp53/tp53_w_1000.fasta"

params.full_reference_genome = "/home/dami/data/ref_genomes/hg_19/hg19.p13.plusMT.full_analysis_set.fa*"

// params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
// if (params.gtf) {
//     Channel.fromPath(params.gtf)
//            .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
//            .into { my_gtf_channel }
// }
/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
indexed_ref = Channel.fromPath(params.full_reference_genome)
// println indexed_ref.view()
input_reads = Channel.fromPath(params.input_reads)
// input_reads_count = input_reads.count().view( x -> "found $x read files")
// input_reads.first().view()

log.info """
    ===================================================
    Cyclomics/CycloSeq : Cyclomics consensus pipeline
    ===================================================
        input_reads : $params.input_reads
        backbone_locus : $params.backbone_locus
        backbone_fasta : $params.backbone_fasta
        backbone_prime_length : $params.backbone_prime_length
        insert_target : $params.insert_target
        insert_target_genome : $params.insert_target_genome
        full_reference_genome : $params.full_reference_genome
        genome               : $params.genome
        
"""

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

//
// WORKFLOW: Run main analysis pipeline
//

include {BwaIndex;
        BwaMemSorted
} from "./modules/bwa.nf"

include {
    Extract5PrimeFasta;
    Extract3PrimeFasta;
} from "./modules/primes.nf"

include {
    Tidehunter;
    TidehunterFullLength;
    TidehunterAggregate;
    Tidehunter53;
    TideHunterTrimmmerFasta;
    TideHunterTrimmmerPrimer;
} from "./modules/tidehunter.nf"


workflow {
    Extract3PrimeFasta(params.backbone_fasta, params.backbone_prime_length)
    Extract5PrimeFasta(params.backbone_fasta, params.backbone_prime_length)
    
    Tidehunter(input_reads)
    TidehunterFullLength(input_reads,
        Extract3PrimeFasta.out,
        Extract5PrimeFasta.out)

    TidehunterAggregate(Tidehunter.out,
    TidehunterFullLength.out)
    TideHunterTrimmmerPrimer(TidehunterAggregate.out)

    // Tidehunter53(
    //     input_reads,
    //     Extract3PrimeFasta.out,
    //     Extract5PrimeFasta.out
    // )
//     Tidehunter(
//         input_reads
//     )

//     // BwaIndex(
//     //     input_reads_count
//     // )
// // align consensus to the expected target
//     BwaMemSorted(
//         Tidehunter.out.consensus,
//         indexed_ref,
//     )
}

/*
========================================================================================
    THE END
========================================================================================
*/
