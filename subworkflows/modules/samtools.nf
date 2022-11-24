#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SamtoolsIndex{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'


    input:
        each path(bam)

    output:
        tuple path(bam), path("*.bai")

    script:
        """
        samtools index $bam
        """
}

process SamtoolsIndexWithID{
    label 'many_cpu_medium'

    input:
        tuple val(x), path(input_bam) 

    output:
       tuple val(x), path(input_bam), path("*.bai")

    script:
        """
        samtools index $input_bam
        """
}

process SamtoolsSort{
    // Given a sam or bam file, make a sorted bam
    // Does break when sam is not `propper` eg no @SQ tags
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path(input_sam)

    output:
        path "${input_sam.SimpleName}_sorted.bam"

    script:
        """
        samtools sort $input_sam -O bam -o "${input_sam.SimpleName}_sorted.bam"
        """
}

process SamtoolsDepth{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(input_bam), path(bai)

    output:
        path "${input_bam.SimpleName}.bed"

    script:
        """
        samtools depth $input_bam > ${input_bam.SimpleName}.bed
        """
}

process SamtoolsDepthToTSV{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // the slim version of buster is missing ps, which is needed for nextflow
    label 'many_med_cpu_huge_mem'

    input:
        path(input_bed)

    output:
        path "${bam_in.simpleName}.depth.tsv"

    script:
        """
        cat $input_bed | sort -r -n -k3 | awk '!x[\$1]++' > ${input_bed.SimpleName}.depth.tsv
        """
}


process SamToBam{
    // Sort, convert and index 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(sam_file)

    output:
        tuple val(X), path("${X}.bam"), path("${X}.bam.bai")

    script:
        """
        samtools view -Sb $sam_file | samtools sort /dev/stdin -o ${X}.bam
        samtools index ${X}.bam
        """
}

process SamtoolsMerge{
    //  samtools merge – merges multiple sorted files into a single file 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path(input_sam)
        val output_base

    output:
        path "${output_base}.bam"

    script:
        """
        samtools merge -O bam ${output_base}.bam $input_sam
        """
}

process SamtoolsQuickcheck{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        stdout

    script:
        """
        samtools quickcheck $bam_in && echo 'Samtools quickcheck ok' || echo 'Samtools quickcheck fail!'
        """
}

process SamtoolsFlagstats{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)
    
    output:
        path("${bam_in.SimpleName}.flagstats_metadata.json")

    script:
        // TODO: get all parameters available
        """
        samtools flagstat -O json $bam_in > ${bam_in.SimpleName}.flagstats.json 
        jq '.["QC-passed reads"] | {additional_info: {"Reference_aligned_with_backbone":."primary mapped"}}' ${bam_in.SimpleName}.flagstats.json > ${bam_in.SimpleName}.flagstats_metadata.json
        """
}
process SamtoolsIdxStats{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)
    
    output:
        path("${bam_in.SimpleName}.idxstats_metadata.json")

    script:
        """
        TOTALREFREADS=\$(samtools idxstats $bam_in | grep -v "^BB" | grep -v "^*" | awk -F' ' '{sum+=\$3;} END{print sum;}')
        echo \$TOTALREFREADS
        TOTALREFREADS=\$TOTALREFREADS jq -n '{additional_info:{"total_reference_mapping_reads":env.TOTALREFREADS,}}' > ${bam_in.SimpleName}.idxstats_metadata.json

        """
}

process SamtoolsFlagstatsMapPercentage{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam), path(bai)
    
    output:
        tuple val(X), stdout

    script:
        // samtools flagstat $bam | grep % | cut -d '(' -f 2 | cut -d '%' -f1
        """
        set -e -o pipefail
        samtools flagstat $bam | grep % | cut -d '(' -f 2 | cut -d '%' -f1 | head -n 1
        """
}

process SamtoolsMergeTuple{
    //  merge n number of bams into one
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path("${X}.merged.bam"), path("${X}.merged.bam.bai")
    
    script:
    // contains bug, remove in future release
    println('Warning Depreciated')
    """
    samtools merge -O bam ${X}.merged.bam \$(find . -name '*.bam')
    samtools index ${X}.merged.bam
    """
}

process SamtoolsMergeBams{
    //  merge n number of bams into one
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        val(X)
        file(bam_in)

    output:
        tuple val(X), path("${X}.merged.bam"), path("${X}.merged.bam.bai")
    
    script:
    """
    ls
    samtools merge -O bam ${X}.merged.bam \$(find . -name '*.bam')
    samtools index ${X}.merged.bam
    """
}

process RemoveUnmappedReads{
    // Sort, convert and index 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.mapped.bam"), path("${X}.mapped.bam.bai")

    script:
        """
        samtools view -b -F 4 $bam_in > ${X}.mapped.bam
        samtools index ${X}.mapped.bam
        """
}

process PrimaryMappedFilter{
    // Sort, convert and index 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.primary_mapped.bam"), path("${X}.primary_mapped.bam.bai")

    script:
        """
        samtools view -b -F 256 $bam_in > ${X}.primary_mapped.bam
        samtools index ${X}.primary_mapped.bam
        """
}

process MapqAndNMFilter{
    // Sort, convert and index 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    
    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.NM_50_mapq_20.bam"), path("${X}.NM_50_mapq_20.bam.bai")

    script:
        """
        samtools view -b -o ${X}.NM_50_mapq_20.bam $bam_in --input-fmt-option 'filter=[NM]<50 && mapq >20'
        samtools index ${X}.NM_50_mapq_20.bam
        """
}

process FindRegionOfInterest{
    //  publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    
    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}_roi.bed")

    script:
        """
        samtools depth $bam_in | awk '\$3>100' | awk '{print \$1"\t"\$2"\t"\$2 + 1}' | bedtools merge -d 25 -i /dev/stdin > ${X}_roi.bed
        """
}

process MPileup{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    
    input:
        tuple val(X), path(bam_in), path(bai_in), path(reference)
        tuple val(X), path(bed)

    output:
        //  use simplename as X and the bed are not unique
        path ("${bam_in.simpleName}.pileup")

    script:
        """
        samtools mpileup $bam_in -l $bed -f $reference --reverse-del --output-QNAME --output-MQ --max-depth 0 > ${bam_in.simpleName}.pileup
        """
}

process BamTagFilter{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    
    input:
        tuple val(X), path(bam_in), path(bai_in)
        val(tag)
        val(minimum_repeats)

    output:
        tuple val(X), path("${X}.${tag}_gt_${minimum_repeats}.bam"), path("${X}.${tag}_gt_${minimum_repeats}.bam.bai")

    script:
        """
        samtools view -b -o ${X}.${tag}_gt_${minimum_repeats}.bam $bam_in --input-fmt-option 'filter=[${tag}]>${minimum_repeats}'
        samtools index ${X}.${tag}_gt_${minimum_repeats}.bam
        """
}
