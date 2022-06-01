


process AddDepthToJson{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'python:3.8'
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

process AnnotateBam{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'quay.io/biocontainers/pysam:0.19.0--py39h5030a8b_0'
    label 'many_cpu_intensive'

    input:
        tuple val(X), path(bam), path(bai)
        path(sequencing_summary)

    output:
        tuple val(X), path("${X}.tag_annotated.sort.bam"), path("${X}.tag_annotated.sort.bam.bai")

    script:
        """
        annotate_bam.py $sequencing_summary $bam ${X}.tag_annotated
        """
}

process AnnotatePartialBam{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'quay.io/biocontainers/pysam:0.19.0--py39h5030a8b_0'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(bam), path(bai), path(sequencing_summary)
        

    output:
        tuple val(X), path("${X}.tag_annotated.sort.bam"), path("${X}.tag_annotated.sort.bam.bai")

    script:
        """
        annotate_bam.py $sequencing_summary $bam ${X}.tag_annotated
        """
}