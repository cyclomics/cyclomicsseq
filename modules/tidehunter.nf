process Tidehunter{
    // _tide_consensus.fasta in ConCall
    publishDir "$baseDir/data/out/$workflow.runName/tidehunter"
    container 'tidehunter'
    
    input:
        path fasta

    output:
        path "${fasta.SimpleName}.consensus.fasta", emit: consensus

    script:
        """
        TideHunter -t ${task.cpus} $fasta > ${fasta.SimpleName}.consensus.fasta
        """
}

process TidehunterPrimersFixedParams{
    // tide_consensus_full_length in ConCall
    publishDir "$baseDir/data/out/$workflow.runName/tidehunter"
    container 'tidehunter'
    
    input:
        path fasta
        path prime_3
        path prime_5

    output:
        path "${fasta.SimpleName}.consensus.fasta", emit: consensus

    script:
        """
        TideHunter -t ${task.cpus} -5 $prime_5 -3 $prime_3 -p 20 -a 0.60 > ${fasta.SimpleName}.consensus.fasta
        """
}

process Tidehunter53{
    publishDir "$baseDir/data/out/$workflow.runName/tidehunter"
    container 'tidehunter'
    
    input:
        path fasta
        path prime_3
        path prime_5

    output:
        path "${fasta.SimpleName}.consensus.tsv", emit: consensus

    script:
        """
        echo "$params.tidehunter.headerlines" > ${fasta.SimpleName}.consensus.tsv
        TideHunter -f 2 -t ${task.cpus} -5 $prime_5 -3 $prime_3 -p $params.tidehunter.minimum_period -a $params.tidehunter.minimum_match_ratio $fasta >> ${fasta.SimpleName}.consensus.tsv
        """
}

process TideHunterTrimmmerFasta {
    // Plug tidehunter fasta into this
    publishDir "$baseDir/data/out/$workflow.runName/tidehunter/trimmed"

    input:
        path original_fasta
    
    output:
        path "${original_fasta.SimpleName}.full_length.fasta", emit: fasta_full_length 

    script:
        """
        sed 's/,.*//' $original_fasta > ${original_fasta.SimpleName}.full_length.fasta
        """
}   

process TideHunterTrimmmerPrimer {
    // plug tidehunter primer into this
    publishDir "$baseDir/data/out/$workflow.runName/tidehunter/trimmed"

    input:
        path consensus_fasta
    
    output:
        path "${consensus_fasta.SimpleName}_consensus_full_length.fasta", emit: fasta_all 
        path "${consensus_fasta.SimpleName}metadata.txt", emit: metadata_all 

    script:
        // # only export metadata from flexible parameter since it is the same if the pattern is found.
        """
        sed -n -e 's/^>//p' $consensus_fasta > ${consensus_fasta.SimpleName}metadata.txt
        sed 's/,.*//' $consensus_fasta > ${consensus_fasta.SimpleName}_consensus_full_length.fasta
        """
}   

