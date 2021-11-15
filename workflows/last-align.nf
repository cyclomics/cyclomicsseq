#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads            = "/home/dami/projects/ENZA004_spin_data_big_insert/concal_results/9f7eb81b-f53c/05_aggregated/9f7eb81b-f53c_tide_consensus_trimmed.fasta"
params.reference_genome = "//home/dami/projects/ENZA004_spin_data_big_insert/reference/Viroflay_ref_F7-R8.fa"
params.output_dir = "data/out/$workflow.runName"

// can be trained or not-trained, and used with fastqs or fasta files
params.method = "trained"  // trained or non-trained
params.input_type = "fastq"  // fastq or fasta

include {
    LastalAlignFasta
    LastalAlignFastq
    LastalAlignTrainedFasta
    LastalAlignTrainedFastq
} from "../subworkflows/align.nf"


reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
ref_ch = Channel.fromPath(params.reference_genome, checkIfExists: true)

workflow {
    aligned = Channel.empty()
    if (params.method == "trained"){
        if (params.input_type == "fasta"){
            aligned = LastalAlignTrainedFasta(
                reads_ch,
                ref_ch
            )
        }
        else if (params.input_type == "fastq"){
            aligned = LastalAlignTrainedFastq(
                reads_ch,
                ref_ch
            )
        }
    }
    if (params.method == "non-trained"){
        if (params.input_type == "fasta"){
            aligned = LastalAlignFasta(
                reads_ch,
                ref_ch
            )
        }
        else if (params.input_type == "fastq"){
            aligned = LastalAlignFastq(
                reads_ch,
                ref_ch
            )
        }
    }
    emit:
        aligned
     
}