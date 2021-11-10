
nextflow.enable.dsl=2

process PYCOQC {
    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"

    input:
        path summary

    output:
        path "*.html" , emit: html
        path "*.json" , emit: json

    script:
    """
    pycoQC \\
        -f $summary \\
        -o pycoqc.html \\
        -j pycoqc.json
    """
}
