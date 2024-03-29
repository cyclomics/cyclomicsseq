process AddDepthToJson{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)
        path(depth_json)

    output:
        tuple val(X), path("${X}.depth.json")

    script:
        """  
        add_depth_info_json.py --global_json $tidehuntertable --depth_json $depth_json --output ${X}.depth.json
        """
}

process AnnotateBamXTags{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // publishDir "${params.output_dir}/consensus_aligned_tagged", mode: 'copy'
    label 'few_very_memory_intensive'

    input:
        tuple val(X), path(bam), path(bai)
        path(sequencing_summary)

    output:
        tuple val(X), path("${X}.tagged.bam"), path("${X}.tagged.bam.bai")

    script:
        """
        annotate_bam_x.py $sequencing_summary $bam ${X}.tagged.bam
        """
}

process AnnotateBamYTags{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(bam), path(json)

    output:
        tuple val(X), path("${X}.annotated.bam"), path("${X}.annotated.bam.bai")

    script:
        """
        samtools index $bam
        annotate_bam_y.py $json $bam ${X}.annotated.bam
        """
}

process CollectClassificationTypes{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        path(metadata_json)

    output:
        path("classification_count.txt")
    
    script:
        """
        gather_readtypes.py "*.metadata.json" classification_count.txt
        """
}

process FindVariants{
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'max_performance'

    input:
        path(reference_genome)
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)

    output:
        tuple path("${bam.simpleName}.snp.vcf"), path("${bam.simpleName}.indel.vcf")
    
    script:
        // We sleep and access the reference genome, since in some rare cases the file needs accessing to 
        // not cause issues in the python code. 
        """
        sleep 1
        ls
        head $reference_genome
        determine_vaf.py $reference_genome $validation_bed $bam ${bam.simpleName}.snp.vcf ${bam.simpleName}.indel.vcf --threads ${task.cpus} 
        """
}

process FilterValidateVariants{
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(snp_vcf), path(indel_vcf), val(X), path(perbase_table)

    output:
        tuple path("${snp_vcf.simpleName}.filtered.snp.vcf"), path("${indel_vcf.simpleName}.filtered.indel.vcf")
    
    script:
        """
        vcf_filter.py -i $snp_vcf -o ${snp_vcf.simpleName}.filtered.snp.vcf \
        --perbase_table $perbase_table \
        --dynamic_vaf_params $params.dynamic_vaf_params_file \
        --min_ao $params.snp_filters.min_ao \
        --min_dpq $params.snp_filters.min_dpq \
        --min_dpq_n $params.snp_filters.min_dpq_n \
        --min_dpq_ratio $params.snp_filters.min_dpq_ratio \
        --max_sap 0 \
        --min_rel_ratio $params.snp_filters.min_rel_ratio \
        --min_abq $params.snp_filters.min_abq

        vcf_filter.py -i $indel_vcf -o ${indel_vcf.simpleName}.filtered.indel.vcf \
        --perbase_table $perbase_table \
        --dynamic_vaf_params $params.dynamic_vaf_params_file \
        --min_ao $params.indel_filters.min_ao \
        --min_dpq $params.indel_filters.min_dpq \
        --min_dpq_n $params.indel_filters.min_dpq_n \
        --min_dpq_ratio $params.indel_filters.min_dpq_ratio \
        --max_sap 0 \
        --min_rel_ratio $params.indel_filters.min_rel_ratio \
        --min_abq $params.indel_filters.min_abq   
        """
}

process MergeNoisyVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(noisy_snp_vcf), path(noisy_indel_vcf)

    output:
        path("${noisy_snp_vcf.simpleName}.vcf")
    
    script:
        """
	    mkdir tmpdir

        bcftools sort $noisy_snp_vcf -o ${noisy_snp_vcf.simpleName}.sorted.snp.vcf --temp-dir tmpdir
        bgzip ${noisy_snp_vcf.simpleName}.sorted.snp.vcf
        tabix ${noisy_snp_vcf.simpleName}.sorted.snp.vcf.gz

        bcftools sort $noisy_indel_vcf -o ${noisy_indel_vcf.simpleName}.sorted.indel.vcf --temp-dir tmpdir
        bgzip ${noisy_indel_vcf.simpleName}.sorted.indel.vcf
        tabix ${noisy_indel_vcf.simpleName}.sorted.indel.vcf.gz

        bcftools concat -a ${noisy_snp_vcf.simpleName}.sorted.snp.vcf.gz ${noisy_indel_vcf.simpleName}.sorted.indel.vcf.gz \
        -O v -o ${noisy_snp_vcf.simpleName}.tmp.vcf
        bcftools sort ${noisy_snp_vcf.simpleName}.tmp.vcf -o ${noisy_snp_vcf.simpleName}.vcf --temp-dir tmpdir
        rm ${noisy_snp_vcf.simpleName}.tmp.vcf

	    rm -r tmpdir
        """
}

process MergeFilteredVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(filtered_snp_vcf), path(filtered_indel_vcf)

    output:
        path("${filtered_snp_vcf.simpleName}_filtered.vcf")
    
    script:
        """
	    mkdir tmpdir

        bcftools sort $filtered_snp_vcf -o ${filtered_snp_vcf.simpleName}_filtered.sorted.snp.vcf --temp-dir tmpdir
        bgzip ${filtered_snp_vcf.simpleName}_filtered.sorted.snp.vcf
        tabix ${filtered_snp_vcf.simpleName}_filtered.sorted.snp.vcf.gz

        bcftools sort $filtered_indel_vcf -o ${filtered_indel_vcf.simpleName}_filtered.sorted.indel.vcf --temp-dir tmpdir
        bgzip ${filtered_indel_vcf.simpleName}_filtered.sorted.indel.vcf
        tabix ${filtered_indel_vcf.simpleName}_filtered.sorted.indel.vcf.gz

        bcftools concat -a ${filtered_snp_vcf.simpleName}_filtered.sorted.snp.vcf.gz ${filtered_indel_vcf.simpleName}_filtered.sorted.indel.vcf.gz \
        -O v -o ${filtered_snp_vcf.simpleName}_filtered.tmp.vcf
        bcftools sort ${filtered_snp_vcf.simpleName}_filtered.tmp.vcf -o ${filtered_snp_vcf.simpleName}_filtered.vcf --temp-dir tmpdir
        rm ${filtered_snp_vcf.simpleName}_filtered.tmp.vcf

	    rm -r tmpdir
        """
}

process AnnotateVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        path(variant_vcf)


    output:
        path("${variant_vcf.simpleName}_annotated.vcf")
    
    script:
        """
        annotate_vcf.py $variant_vcf ${variant_vcf.simpleName}_annotated.vcf
        """
}

process PlotFastqsQUalAndLength{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
        path(fastq)
        val grep_pattern
        val plot_file_prefix
        val tab_name

    output:
        path("${plot_file_prefix}_histograms.json")
    
    script:
        """
        plot_fastq_histograms.py $grep_pattern ${plot_file_prefix}_histograms.html $tab_name $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${plot_file_prefix}_histograms.json
        """
}

process PlotReadStructure{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_cpu_medium'

    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        tuple val(X), path(bam), path(bai)

    output:
        path("${bam.simpleName}_read_structure.json")

    script:
        // This takes a lot of RAM when the sequencing summary is big!
        """
        samtools sort -n -o tmp_readname_sorted_${bam.simpleName}.bam ${bam}
        plot_read_structure_donut.py tmp_readname_sorted_${bam.simpleName}.bam \
        ${bam.simpleName}_read_structure.json \
        $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${bam.simpleName}_read_structure.json
        """
}

process PlotVcf{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        path(vcf)

    output:
        path("${vcf.simpleName}.json")
    
    script:
        // This takes a lot of RAM when the sequencing summary is big!
        """
        plot_vcf.py $vcf ${vcf.simpleName}.html $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${vcf.simpleName}.json
        """
}

process PasteVariantTable{
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
        path(vcf_file)

    output:
        path("${vcf_file.simpleName}_table.json")
    
    script:
        """
        write_variants_table.py $vcf_file ${vcf_file.simpleName}_table.json 'Variant table' $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${vcf_file.simpleName}_table.json
        """
}

process PlotQScores{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
        tuple val(X), path(split_pileup)
        tuple val(Y), path(consensus_pileup)

    output:
        path("${consensus_pileup.simpleName}.json")
    
    script:
        """
        plot_bam_accuracy.py $split_pileup $consensus_pileup ${consensus_pileup.simpleName}.html $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${consensus_pileup.simpleName}.json
        """
}

process PlotMetadataStats{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_huge_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        path(jsons)

    output:
        path("metadata_plots.json")
    
    script:
        """
        plot_metadata.py . metadata_plots.html $params.priority_limit --subsample_size $params.metadata.subsample_size \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt metadata_plots.json
        """
}

process PlotReport{
    publishDir "${params.output_dir}", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        path(jsons)

    output:
        path("report.html")
    
    script:
        """
        generate_report.py '${params}' $workflow.manifest.version $params.priority_limit
        """
}