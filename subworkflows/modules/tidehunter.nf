#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Tidehunter{
    // _tide_consensus.fasta in ConCall
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'quay.io/biocontainers/tidehunter:1.5.3--h2e03b76_0'
    label 'many_med_cpu_huge_mem'

    input:
        path fasta

    output:
        path "${fasta.SimpleName}_tide_consensus.fasta", emit: consensus

    script:
        """
        TideHunter -t ${task.cpus} $fasta > ${fasta.SimpleName}_tide_consensus.fasta
        """
}

process TidehunterLongest{
    // _tide_consensus.fasta in ConCall
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'quay.io/biocontainers/tidehunter:1.5.3--h2e03b76_0'
    label 'many_med_cpu_huge_mem'

    input:
        path fasta

    output:
        path "${fasta.SimpleName}_tide_consensus.tsv", emit: consensus

    script:
        """
        TideHunter -f 2 -t ${task.cpus} -k ${params.tidehunter.kmer_length} --min-copy ${params.tidehunter.min_copy} --min-period ${params.tidehunter.min_period} --longest --thread ${task.cpus} $fasta > ${fasta.SimpleName}_tide_consensus.tsv
        """
}

process TideHunterTableToFasta{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path tidehuntertable
    
    output:
        path "${tidehuntertable.SimpleName}.fasta", emit: consensus

    script:
        """
        awk '{print ">"\$1"_"\$3"_2_3_4", "\\n", \$11}' $tidehuntertable > ${tidehuntertable.SimpleName}.fasta
        """
}

process TideHunterTrimmmer {
    // Plug tidehunter fasta into this
    // removes everyting after the first comma, since this messes with the bam spec
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path original_fasta
    
    output:
        path "${original_fasta.SimpleName}.trimmed.fasta", emit: trimmed_fasta 
        path "${original_fasta.SimpleName}.original_fasta_names.txt", emit: original_fasta_names

    script:
        """
        cut -d ' ' -f 1 $original_fasta > ${original_fasta.SimpleName}.trimmed.fasta
        grep '>' $original_fasta > ${original_fasta.SimpleName}.original_fasta_names.txt
        """
}   


process Tidehunter53QualTable{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'quay.io/biocontainers/tidehunter:1.5.3--h2e03b76_0'
    label 'many_cpu_medium'


    input:
        tuple path(fasta),  path(prime_3_fasta),  path(prime_5_fasta)

    output:
        tuple( val("${fasta.baseName}"), path("*consensus.tsv"))

    script:
        """
        echo "$params.tidehunter.headerlinesQual" > ${fasta.SimpleName}.consensus.tsv
        TideHunter \
        --out-fmt 4  \
        --kmer-length 16 \
        --longest  \
        --thread ${task.cpus} \
        --five-prime $prime_5_fasta \
        --three-prime $prime_3_fasta \
        --min-period $params.tidehunter.minimum_period \
        --min-len $params.tidehunter.minimum_length \
        --min-copy $params.tidehunter.minimum_copy \
        -a $params.tidehunter.minimum_match_ratio \
        $fasta >> ${fasta.SimpleName}.consensus.tsv
        """
}
process TideHunterFilterTableStartpos{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple( val(X), path(tidehuntertable))
    
    output:
        tuple( val(X), path("${tidehuntertable.baseName}.startposfilter.tsv"))

    script:
        """
            awk 'BEGIN {getline; print }{ if (\$5 < 100) { print }}' $tidehuntertable > ${tidehuntertable.baseName}.startposfilter.tsv
        """
}

process TideHunterQualTableToFastq{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)

    
    output:
        tuple val(X), path("${tidehuntertable.SimpleName}.fastq")

    script:
        //  skip first line (header), go on to convert all lines to fastq entries
        """
            awk 'BEGIN { getline }{print "@"\$1"\\n"\$11"\\n+\\n"\$12}' $tidehuntertable > ${tidehuntertable.SimpleName}.fastq
        """
}

process TideHunterQualTableToJson{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)

    output:
        tuple val(X), path("${X}.json")

    script:
        """
        jq --raw-input --slurp \
        'split("\n") | 
        map(split("\t")) | 
        .[0:-1] | 
        map( { "id": .[0], 
            "raw_length":.[3],  
            "baseunit_copies": .[2], 
            "baseunit_start_idx":.[4],  
            "baseunit_end_idx":.[5],  
            "baseunit_length":.[6],
            "baseunit_certainty":.[7],
            "baseunit_orientation":.[8],
            "baseunit_idx":.[9],  
        })' \
        $tidehuntertable > ${X}.json
        """
}

process TideHunterQualJsonMerge{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        val(X)
        path(tidehuntertable)

    output:
        tuple val(new_X), path("${new_X}.json")

    script:
        new_X = X.split('_')[0]
        // remove all items containing the column names
        """
        jq  'flatten | .[] | select(.id=="readName" | not)' *.json | jq --slurp '.' > ${new_X}.json
        """
}
