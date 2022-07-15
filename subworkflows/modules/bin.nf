


process AddDepthToJson{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
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
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(metadata_json)

    output:
        path("classification_count.txt")
    
    script:
        """
        gather_readtypes.py "*.metadata.json" classification_count.txt
        """
}

process VariantValidate{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)


    output:
        path("${bam.simpleName}_validated.vcf")
    
    script:
        """
        determine_vaf.py $validation_bed $bam ${bam.simpleName}_validated.vcf
        """
}

process FilterVariants {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(vcf_file)

    output:
        path("${vcf_file.simpleName}_filtered.vcf")
    
    script:
        """
        vcf_filter.py $vcf_file ${vcf_file.simpleName}_filtered.vcf
        """
}

process PlotFastqsQUalAndLength{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(fastq)
        val grep_pattern
        val plot_file_prefix

    output:
        path("${plot_file_prefix}_histograms.html")
    
    script:
        """
        plot_fastq_histograms.py $grep_pattern ${plot_file_prefix}_histograms.html
        """
}

process PlotReadStructure{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        path("${bam.simpleName}_aligned_segments.html")
        path("${bam.simpleName}_read_structure.html")

    
    script:
        """
        plot_read_structure_donut.py $bam ${bam.simpleName}_aligned_segments.html ${bam.simpleName}_read_structure.html
        """
}

process PlotVcf{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(vcf)

    output:
        path("${vcf.simpleName}.html")
    
    script:
        """
        plot_vcf.py $vcf ${vcf.simpleName}.html
        """
}

process PlotQScores{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        tuple val(X), path(split_pileup)
        tuple val(Y), path(consensus_pileup)

    output:
        path("${consensus_pileup.simpleName}.html")
    
    script:
        """
        
        plot_bam_accuracy.py $split_pileup $consensus_pileup ${consensus_pileup.simpleName}.html
        """
}

process PlotMetadataStats{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        val(X)
        path(jsons)

    output:
        path("metadata_plots.html")
    
    script:
        """
        plot_metadata.py . metadata_plots.html
        """
}
