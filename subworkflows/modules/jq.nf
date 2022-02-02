
process JqAddDepthToJson{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'stedolan/jq'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)
        path(depth_json)

    output:
        tuple val(X), path("${X}.depth.json")

    script:
        """
        jq --argjson text "\$(jq -c '' $depth_json )" \
        '.depth=\$text' global.json > ${X}.depth.json
        """
}
