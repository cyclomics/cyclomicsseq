#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MedakaSmolecule {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'ontresearch/medaka:v1.6.0'
    label 'many_low_cpu_tiny_mem'

    // errorStrategy 'ignore' 

    input:
        path(input_tar)

    output:
        tuple val("${input_tar.simpleName}"), path("${input_tar.simpleName}.consensus.fastq")

    script:
        //  We hide the for loop in a script so that we can catch exit code 1 cases if medaka fails to make consensus
        """
        tar xvf $input_tar
        mkdir /medaka_results
        bash medakasmolecule.sh ${params.medaka.depth} ${params.medaka.model} ${params.medaka.chunk_len} ${params.medaka.chunk_ovlp} ${params.medaka.method} ${params.medaka.length} ${task.cpus}  ${params.medaka.batch_size} ${input_tar.simpleName}.consensus.fastq
        echo "done with medaka loop"
        cat /medaka_results/*/consensus.fastq > ${input_tar.simpleName}.consensus.fastq
        """
        // for i in \$(ls cycas_results); do
        //     echo \$i
        //     medaka smolecule --depth 1 --model r104_e81_sup_g5015 --chunk_len 100 --chunk_ovlp 50 --method spoa --length 50 --threads 1 --qualities --batch_size 1000 /medaka_results/ cycas_results/\$i
        //     mv /medaka_results/consensus.fastq \${i}.consensus.fastq
        // done
}
