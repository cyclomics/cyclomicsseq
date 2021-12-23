
nextflow.enable.dsl=2

process Tidehunter{
    // _tide_consensus.fasta in ConCall
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'tidehunter'
    
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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'tidehunter:1.5.1'
    
    cpus = 1

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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'tidehunter:1.5.3'
    
    cpus = 2

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
        --min-len $params.tidehunter.minimum_period \
        --min-copy $params.tidehunter.minimum_copy \
        -a $params.tidehunter.minimum_match_ratio \
        $fasta >> ${fasta.SimpleName}.consensus.tsv
        """
}
process TideHunterFilterTableStartpos{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

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
