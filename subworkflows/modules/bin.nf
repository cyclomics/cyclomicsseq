


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
        tuple val(X), path("${X}.annotated.bam"), path("${X}.annotated.bam.bai")

    script:
        """
        annotate_bam_x.py $sequencing_summary $bam ${X}.annotated.bam
        """
}

process AnnotateBamYTags{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(bam), path(bai), path(json)
        

    output:
        tuple val(X), path("${X}.y_tagged.sort.bam"), path("${X}.y_tagged.sort.bam.bai")

    script:
        """
        annotate_bam_y.py $json $bam ${X}.y_tagged
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
        path(validation_bed)


    output:
        path("${bam.simpleName}_validated.vcf")
    
    script:
        """
        determine_vaf.py $validation_bed $bam ${bam.simpleName}_validated.vcf
        """
}
