


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
    publishDir "${params.output_dir}/consensus_aligned", mode: 'copy'

    label 'many_cpu_intensive'

    input:
        tuple val(X), path(bam), path(bai)
        path(sequencing_summary)

    output:
        tuple val(X), path("${X}.taged.bam"), path("${X}.taged.bam.bai")

    script:
        """
        annotate_bam_x.py $sequencing_summary $bam ${X}.taged.bam
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

    input:
        path(metadata_json)

    output:
        path("classification_count.txt")
    
    script:
        """
        gather_readtypes.py "*.metadata.json" classification_count.txt
        """
}

process FindSNPs{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)

    output:
        path("${bam.simpleName}.snp.vcf")
    
    script:
        """
        detect_snp.py $validation_bed $bam ${bam.simpleName}.snp.vcf
        """
}

process FilterSNPs {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        path(vcf_file)

    output:
        tuple path("${vcf_file.simpleName}_filtered.snp.vcf.gz"), path("${vcf_file.simpleName}_filtered.snp.vcf.gz.tbi")
    
    script:
        """
        vcf_filter.py $vcf_file ${vcf_file.simpleName}_filtered.snp.tmp.vcf
        bcftools sort ${vcf_file.simpleName}_filtered.snp.tmp.vcf -o ${vcf_file.simpleName}_filtered.snp.vcf
        rm ${vcf_file.simpleName}_filtered.snp.tmp.vcf
        bgzip ${vcf_file.simpleName}_filtered.snp.vcf
        tabix ${vcf_file.simpleName}_filtered.snp.vcf.gz
        """
}

process FindIndels{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        path(reference_genome)
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)

    output:
        tuple path("${bam.simpleName}.indel.vcf.gz"), path("${bam.simpleName}.indel.vcf.gz.tbi")
    
    script:
        """
        detect_indel.py $reference_genome $validation_bed $bam ${bam.simpleName}.indel.tmp.vcf
        bcftools sort ${bam.simpleName}.indel.tmp.vcf -o ${bam.simpleName}.indel.vcf
        rm ${bam.simpleName}.indel.tmp.vcf
        bgzip ${bam.simpleName}.indel.vcf
        tabix ${bam.simpleName}.indel.vcf.gz
        """
}

process MergeVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        tuple path(snp_vcf), path(snp_tbi), path(indel_vcf), path(indel_tbi)

    output:
        path("${indel_vcf.simpleName}.merged.vcf")
    
    script:
        """
        bcftools concat -a $snp_vcf $snp_vcf -O v -o ${indel_vcf.simpleName}.merged.tmp.vcf
        bcftools sort ${indel_vcf.simpleName}.merged.tmp.vcf -o ${indel_vcf.simpleName}.merged.vcf
        rm ${indel_vcf.simpleName}.merged.tmp.vcf
        """
}

process AnnotateVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'

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

    input:
        path(fastq)
        val grep_pattern
        val plot_file_prefix
        val tab_name

    output:
        tuple path("${plot_file_prefix}_histograms.html"), path("${plot_file_prefix}_histograms.json")
    
    script:
        """
        plot_fastq_histograms.py $grep_pattern ${plot_file_prefix}_histograms.html $tab_name
        """
}

process PlotReadStructure{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple path("${bam.simpleName}_aligned_segments.html"), path("${bam.simpleName}_read_structure.html"), path("${bam.simpleName}_read_structure.json")

    script:
        """
        plot_read_structure_donut.py $bam ${bam.simpleName}_aligned_segments.html ${bam.simpleName}_read_structure.html
        """
}

process PlotVcf{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(vcf)

    output:
        tuple path("${vcf.simpleName}.html"), path("${vcf.simpleName}.json")
    
    script:
        """
        plot_vcf.py $vcf ${vcf.simpleName}.html
        """
}

process PlotQScores{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        tuple val(X), path(split_pileup)
        tuple val(Y), path(consensus_pileup)

    output:
        tuple path("${consensus_pileup.simpleName}.html"), path("${consensus_pileup.simpleName}.csv"), path("${consensus_pileup.simpleName}.json")
    
    script:
        """
        plot_bam_accuracy.py $split_pileup $consensus_pileup ${consensus_pileup.simpleName}.html
        """
}

process PlotMetadataStats{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        tuple val(x), path(jsons)

    output:
        tuple path("metadata_plots.html"), path("metadata_plots.json")
    
    script:
        """
        plot_metadata.py . metadata_plots.html
        """
}
process PlotReport{
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(jsons)

    output:
        path("report.html")
    
    script:
        """
        ls
        generate_report.py '${params}' $workflow.manifest.version
        """

}
