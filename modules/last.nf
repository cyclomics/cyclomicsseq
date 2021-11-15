#!/usr/bin/env nextflow

/*
These processes can be used to run Lastal on a docker container called lastal
*/

nextflow.enable.dsl=2

process LastCreateDB{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'lastal'
    
    input:
        path reference_genome

    output:
        path "db.*"

    script:
        // -uNEAR // old params??
        """
        lastdb db $reference_genome
       """
}

process LastTrainModel{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'lastal'
    
    input:
        path read_fq
        path last_db

    output:
        path "train.out"

    script:
        // Add -Q1 to train on fastq files
        """ 
        /lastal_921/scripts/last-train db $read_fq > train.out
        """
}


process LastSplit{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'lastal'
    
    input:
        path lastal_file

    output:
        path "*.maf"

    script:
        """
        last-split $lastal_file > "${lastal_file.simpleName}.maf"
        """
}

process LastAlignTrained{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'lastal'
    
    input:
        path read_fq
        path last_db
        path trainfile


    output:
        path "*.lastal"

    script:
        // lastal -p $trainfile $last_db $read_fq > "${read_fq.simpleName}.lastal"
        """
        lastal -p $trainfile db $read_fq > "${read_fq.simpleName}.lastal"
        """
}

process LastAlign{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'lastal'

    input:
        path reads
        path db

    output:
        path "${reads.simpleName}.maf"
    script:
        // lastal -p $trainfile $last_db $reads > "${reads.simpleName}.lastal"
        """
        lastal db $reads > ${reads.simpleName}.maf
        """
}

process Maf2sam {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    container "lastal"

    input:
        path(maf)
    
    output:
        path "${maf.SimpleName}.sam"

    script:
    """
    /lastal_921/scripts/maf-convert sam $maf > ${maf.SimpleName}.sam
    """
}

process SamtoolsFixSam{
    // Lastal maf-convert does not add @SQ tags to the sam, this is fixed by using samtools view -T here
    // also sorts and converts to bam
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        path input_sam
        path reference_genome

    output:
        path "${input_sam.SimpleName}_sorted.bam"

    script:
        """
        samtools view -hbS $input_sam -T $reference_genome |samtools sort - -O bam -o "${input_sam.SimpleName}_sorted.bam"
        """

}

