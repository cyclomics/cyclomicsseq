#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MedakaSmolecule {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'ontresearch/medaka:v1.6.0'
    label 'many_low_cpu_tiny_mem'

    input:
        path(input_tar)

    output:
        path("${input_tar.simpleName}.consensus.fastq")

    script:
        """
        tar xvf $input_tar
        mkdir /medaka_results

        for i in \$(ls cycas_results); do
            echo \$i
            medaka smolecule --depth 1 --model r104_e81_sup_g5015 --chunk_len 100 --chunk_ovlp 50 --method spoa --length 50 --threads 1 --qualities --batch_size 1000 /medaka_results/ cycas_results/\$i
            mv /medaka_results/consensus.fastq \${i}.consensus.fastq
        done

        """
            // medaka smolecule \
            //     --depth ${params.medaka.depth} \
            //     --model ${params.medaka.model} \
            //     --chunk_len ${params.medaka.chunk_len} \
            //     --chunk_ovlp ${params.medaka.chunk_ovlp} \
            //     --method ${params.medaka.method} \
            //     --length ${params.medaka.length} \
            //     --threads ${task.cpus} \
            //     --qualities \
            //     --batch_size ${params.medaka.batch_size} \
            //     /medaka_results/ \
            //     cycas_results/$i
            // mv /medaka_results/consensus.fastq $i.consensus.fastq

    // for fq in $(ls cycas_results); do
    //         echo $fq
    //         medaka smolecule \
    //         --depth 1 \
    //         --model r104_e81_sup_g5015 \
    //         --chunk_len 100 \
    //         --chunk_ovlp 50 \
    //         --method spoa \
    //         --length 50 \
    //         --threads 1 \
    //         --qualities \
    //         --batch_size 1000 \ 
    //         /medaka_results/ \
    //         cycas_results/$fq

}
