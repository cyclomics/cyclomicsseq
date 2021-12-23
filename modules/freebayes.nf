process Freebayes {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'nfcore/sarek:2.7.1'
    cpus 1

    input:
        tuple val(X), path(input_bam_file),path(input_bai_file), path(reference)

    output:
        tuple val(X), path("${input_bam_file.SimpleName}.vcf")

    script:
        println("WARN: out chr17 variant calling only!")
        ref = reference.first()
        """
        freebayes \
        --haplotype-length 3 \
        --min-alternate-count 1 \
        -r chr17 \
        -f $ref \
        --vcf ${input_bam_file.SimpleName}.vcf \
        $input_bam_file
        """
}
