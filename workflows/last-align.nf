#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads            = "/home/dami/projects/ENZA004_spin_data_big_insert/concal_results/9f7eb81b-f53c/05_aggregated/9f7eb81b-f53c_tide_consensus_trimmed.fasta"
params.reference_genome = "//home/dami/projects/ENZA004_spin_data_big_insert/reference/Viroflay_ref_F7-R8.fa"
params.output_dir = "data/out/$workflow.runName"

// can be trained or not-trained
params.method = "trained"

include {
    LastalAlign;
    LastalAlignTrained
} from "../subworkflows/align.nf"


reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
ref_ch = Channel.fromPath(params.reference_genome, checkIfExists: true)

workflow {
    
    aligned = LastalAlignTrained(
            reads_ch,
            ref_ch
        )

    emit:
        aligned
    
    
}