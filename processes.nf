#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// REFERENCE AND INPUT PARSING
process Reference_info {
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta

    output:
        path "reference_info.txt"

    script:
        """
        echo "--------- md5 hash info ---------" >> reference_info.txt
        md5sum $fasta >> reference_info.txt
        echo "--------- assembly info ---------" >> reference_info.txt
        seqkit fx2tab --length --name --header-line --seq-hash $fasta >> reference_info.txt
        """
}


process BwaIndex{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_med_cpu_huge_mem'

    input:
        path reference_genome

    output:
        path "${reference_genome}*", emit: bwa_index

    script:
        """
        bwa index $reference_genome
        """
}

process Minimap2Index{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'few_memory_intensive'
    
     input:
        path(reference)
    
    output:
        tuple path(reference), path("${reference}.fai")

    script:
        """
        seqkit faidx $reference
        """
}

process FastaIndex{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'few_memory_intensive'
    
     input:
        path(reference_genome)
    
    output:
        path("${reference_genome.simpleName}.mmi")

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} -d ${reference_genome.simpleName}.mmi $reference_genome
        """
}

process SplitReadFilesOnNumberOfReads {
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(fq)

    output:
        tuple val(sample_id), val(file_id), path("split/${file_id}_*.fastq", arity: '1..*')

    script:
        """
        seqkit split -j ${task.cpus} -s $params.max_fastq_size --by-size-prefix ${file_id}_ -O split $fq
        """
}

process FilterShortReads{
    label 'many_cpu_medium'
    
    input:
        tuple val(sample_id), val(file_id), path(fq)

    output:
        tuple val(sample_id), val("${file_id}_filtered"), path("${file_id}_filtered.fastq")

    script:
        """
        seqkit seq -m ${params.filtering.minimum_raw_length} $fq > ${file_id}_filtered.fastq
        """
}

process MergeFasta {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta1
        path fasta2
    
    output:
        path "${fasta1.simpleName}_${fasta2.simpleName}.fasta"
    
    script:
        """
        cat $fasta1 > ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        cat $fasta2 >> ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        """

}

process SplitReadsOnAdapterSequence {
    // https://github.com/nanoporetech/duplex-tools/blob/master/fillet.md
    label 'many_low_cpu_low_mem'

    input:
        tuple val(sample_id), val(file_id), path(fq)

    output:
        tuple val(sample_id), val(file_id), path("results/${fq.simpleName}_split.fastq.gz")

    script:
        """
        duplex_tools split_on_adapter . results/ Native 
        """
} 



// ALIGNMENT
process SamtoolsMergeBams{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_ids), path(bams_in)

    output:
        tuple val(sample_id), val(sample_id), path("${sample_id}.merged.bam"), path("${sample_id}.merged.bam.bai")
    
    script:
        """
        samtools merge -p -c -O bam ${sample_id}.merged.bam \$(find . -name '*.bam')
        samtools index ${sample_id}.merged.bam
        """
}

process BamTagFilter{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/consensus_aligned", mode: 'copy'
    label 'many_cpu_medium'
    
    input:
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)
        val(tag)
        val(minimum_repeats)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}.${tag}_gt_${minimum_repeats}.bam"), path("${file_id}.${tag}_gt_${minimum_repeats}.bam.bai")

    script:
        """
        samtools view -b -o ${file_id}.${tag}_gt_${minimum_repeats}.bam $bam_in --input-fmt-option 'filter=[${tag}]>=${minimum_repeats}'
        samtools index ${file_id}.${tag}_gt_${minimum_repeats}.bam
        """
}

process SamtoolsIndexWithID{
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(bam)

    output:
        tuple val(sample_id), val(file_id), path(bam), path("*.bai") 

    script:
        """
        samtools index $bam
        """
}

process MapqAndNMFilter{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}.NM_50_mapq_20.bam"), path("${file_id}.NM_50_mapq_20.bam.bai")

    script:
        """
        samtools view -b -o ${file_id}.NM_50_mapq_20.bam $bam_in --input-fmt-option 'filter=[NM]<50 && mapq >20'
        samtools index ${file_id}.NM_50_mapq_20.bam
        """
}

process PrimaryMappedFilter{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}.primary_mapped.bam"), path("${file_id}.primary_mapped.bam.bai")

    script:
        """
        samtools view -b -F 256 $bam_in > ${file_id}.primary_mapped.bam
        samtools index ${file_id}.primary_mapped.bam
        """
}

process Minimap2Align{    
    // Use standard Minimap2 parameters for alignment, also works with .mmi files.
    cpus params.economy_mode == true ? 2 : 7
    
    // *1.8 for T2T gives 13.5 for the first try and 27 for the second
    // Thus fitting within 32GB systems (which have less than 32 in reality)
    // Small refs below 1GB get 4GB per task retry (4,8,12)
    memory = {reference_genome.size() > 1_000_000_000 ? Math.round(reference_genome.size()*1.8 * task.attempt) : "4GB"* task.attempt}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), val(file_id), path(fq)
        path(reference_genome)
    
    output:
        tuple val(sample_id), val(file_id), path("${file_id}.bam")

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fq > tmp.sam 
        samtools sort -o ${file_id}.bam tmp.sam
        rm tmp.sam
        """
}

process Minimap2AlignAdaptiveParameterized{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'minimap_large'

    memory {reference_genome.size() > 31_000_000_000 ? "30GB" : "${reference_genome.size() * (1 + task.attempt)}B"}
    cpus (params.economy_mode == true ? 2 :{reference_genome.size() < 500_000_000 ? 4 : 8 })

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), val(file_id), path(fq)
        path(reference_genome)
    
    output:
        tuple val(sample_id), val(file_id), path("${file_id}.bam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} -m ${params.minimap2parameterized.min_chain_score} -n ${params.minimap2parameterized.min_chain_count} -s ${params.minimap2parameterized.min_peak_aln_score} $reference_genome $fq > tmp.sam 
        samtools sort -o ${file_id}.bam tmp.sam
        rm tmp.sam
        """
}

process AnnotateBamXTags {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // publishDir "${params.output_dir}/consensus_aligned_tagged", mode: 'copy'
    label 'few_very_memory_intensive'

    input:
        tuple val(X), path(bam), path(bai)
        path sequencing_summary

    output:
        tuple val(X), path("${X}.tagged.bam"), path("${X}.tagged.bam.bai")

    script:
        """
        annotate_bam_x.py ${sequencing_summary} ${bam} ${X}.tagged.bam
        """
}

process AnnotateBamYTags {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(sample_id), val(file_id), path(bam), path(json)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}.annotated.bam"), path("${file_id}.annotated.bam.bai")

    script:
        """
        samtools index ${bam}
        annotate_bam_y.py ${json} ${bam} ${file_id}.annotated.bam
        """
}

process BamAlignmentRateFilter {
    publishDir "${params.output_dir}/consensus_aligned", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}.filtered.bam"), path("${file_id}.filtered.bam.bai")

    script:
        """
        filter_alignment_rate.py --input_bam ${bam_in} --output_bam ${file_id}.filtered.bam --amplicons_bed ${params.region_file} --min_align_rate ${params.min_align_rate}
        samtools index ${file_id}.filtered.bam
        """
}

process BwaMemSorted{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

    input:
        each path(fastq)
        file(reference)
        file(reference_indexes)

    output:
        tuple val("${fastq.simpleName}"), val("${fastq.simpleName}"), path("${fastq.simpleName}.bam"), path("${fastq.simpleName}.bam.bai")
        
    script:
        // using grep to find out if ref is fa or fasta, plug in env var to bwa
        // piplefail for better control over failures
        // REMOVED. TODO: check reference file itself
        """
        bwa mem -R "@RG\\tID:${params.bwamem.readgroup}\\tSM:${params.bwamem.sampletag}\\tPL:${params.bwamem.platform}" -M -t ${task.cpus} -c ${params.bwamem.mem_max_genome_occurance} -L ${params.bwamem.softclip_penalty} -M $reference $fastq | \
        samtools sort -@ ${task.cpus} /dev/stdin -o "${fastq.simpleName}.bam"
        samtools index ${fastq.simpleName}.bam
        """
}

process BwaMemContaminants{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_med_cpu_huge_mem'

    input:
        each path(fastq)
        file(reference)
        file(reference_indexes)

    output:
        tuple val("${fastq.simpleName}"), path("${fastq.simpleName}_contaminants.bam"), path("${fastq.simpleName}_contaminants.bam.bai")
        
    script:
        """
        bwa mem -R "@RG\\tID:${params.bwamem.readgroup}\\tSM:${params.bwamem.sampletag}\\tPL:${params.bwamem.platform}" -M -t ${task.cpus} -c ${params.bwamem.mem_max_genome_occurance} -L ${params.bwamem.softclip_penalty} -M $reference $fastq | \
        samtools sort -@ ${task.cpus} /dev/stdin -o "${fastq.simpleName}_contaminants.bam"
        samtools index ${fastq.simpleName}_contaminants.bam
        """
}


// CONSENSUS CALLING
process Cycas{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/consensus", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(bam), path(bai)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}.consensus.fastq"), path("${file_id}.metadata.json")

    script:
        
        """
        mkdir plots
        python $params.cycas_location consensus --bam-file $bam --output ${file_id}.consensus.fastq --metadata-json  ${file_id}.metadata.json
        """
    }

// VARIANT CALLING AND CONTAMINANTION QC
process FindVariants {
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'max_performance'

    input:
        path reference_genome
        tuple val(sample_id), val(file_id), path(bam), path(bai), path(validation_bed)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}_snp.vcf"), path("${file_id}_indel.vcf")

    script:
        // We sleep and access the reference genome, since in some rare cases the file needs accessing to 
        // not cause issues in the python code. 
        """
        sleep 1
        ls
        head ${reference_genome}
        determine_vaf.py ${reference_genome} ${validation_bed} ${bam} ${file_id}_snp.vcf ${file_id}_indel.vcf --threads ${task.cpus} 
        """
}

process FilterValidateVariants {
    //publishDir "${params.output_dir}/variants/filtered_by_mode", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(sample_id), val(file_id), path(vcf), path(perbase_table)
        val mode

    output:
        tuple val(sample_id), val(file_id), path("${file_id}_${mode}_filtered.vcf")

    script:
        if (mode == 'snp') {
            """
            vcf_filter.py -i ${vcf} -o ${file_id}_${mode}_filtered.vcf \
            --perbase_table ${perbase_table} \
            --dynamic_vaf_params ${params.dynamic_vaf_params_file} \
            --min_ao ${params.snp_filters.min_ao} \
            --min_dpq ${params.snp_filters.min_dpq} \
            --min_dpq_n ${params.snp_filters.min_dpq_n} \
            --min_dpq_ratio ${params.snp_filters.min_dpq_ratio} \
            --max_sap 0 \
            --min_rel_ratio ${params.snp_filters.min_rel_ratio} \
            --min_abq ${params.snp_filters.min_abq}
            """
        }
        else if (mode == 'indel') {
            """
            vcf_filter.py -i ${vcf} -o ${file_id}_${mode}_filtered.vcf \
            --perbase_table ${perbase_table} \
            --dynamic_vaf_params ${params.dynamic_vaf_params_file} \
            --min_ao ${params.indel_filters.min_ao} \
            --min_dpq ${params.indel_filters.min_dpq} \
            --min_dpq_n ${params.indel_filters.min_dpq_n} \
            --min_dpq_ratio ${params.indel_filters.min_dpq_ratio} \
            --max_sap 0 \
            --min_rel_ratio ${params.indel_filters.min_rel_ratio} \
            --min_abq ${params.indel_filters.min_abq}   
            """
        }
        else {
            error("Invalid merging mode (VCF type): ${mode}. Valid modes are 'snp' or 'indel'.")
        }
}


process PerbaseBaseDepth {
    publishDir "${params.output_dir}/depth_tables", pattern: "*consensus.tsv", mode: 'copy'
    label 'few_very_memory_intensive'

    input:
        tuple val(sample_id), val(file_id), path(input_bam_file), path(input_bai_file), path(reference)
        tuple val(bed_sample_id), val(bed_file_id), path(bed)
        val(output_name)

    output:
        tuple val(sample_id), val(file_id), path("${sample_id}_${output_name}")

    script:
    """
    samtools faidx $reference
    perbase base-depth -F 256 -D $params.perbase.max_depth -t $task.cpus -z -b $bed --ref-fasta $reference $input_bam_file > ${sample_id}_${output_name}
    """
}



process SortVCF {
    // publishDir "${params.output_dir}/variants/sorted", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(sample_id), val(file_id), path(vcf)

    output:
        tuple val(sample_id), val(file_id), path("${vcf.simpleName}_sorted.vcf.gz"), path("${vcf.simpleName}_sorted.vcf.gz.tbi")

    script:
        """
        bcftools sort ${vcf} -o ${vcf.simpleName}_sorted.vcf --temp-dir .
        bgzip -f ${vcf.simpleName}_sorted.vcf
        tabix -f ${vcf.simpleName}_sorted.vcf.gz
        """
}


process MergeVCF {
    publishDir "${params.output_dir}/variants/noisy", pattern: "*noisy.vcf", mode: 'copy'
    publishDir "${params.output_dir}/variants/filtered", pattern: "*filtered.vcf", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(sample_id), val(file_id), path(snp_vcf), path(snp_tbi), path(indel_vcf), path(indel_tbi)
        val mode

    output:
        tuple val(sample_id), val(file_id), path("${file_id}_${mode}.vcf")

    script:
        if (mode == 'noisy' || mode == 'filtered') {
            """
            bcftools concat -a ${snp_vcf} ${indel_vcf} \
            -O v -o ${file_id}_${mode}.tmp.vcf \
            2> >(tee -a error.txt >&2) && \
            if grep -q "E::" error.txt; then \
            echo "VCF concatenation failed in certain rows.\n"; \
            exit 1; \
            fi

            bcftools sort ${file_id}_${mode}.tmp.vcf -o ${file_id}_${mode}.vcf --temp-dir .
            rm ${file_id}_${mode}.tmp.vcf
            """
        } else {
            error("Invalid merging mode (VCF type): ${mode}. Valid modes are 'noisy' or 'filtered'.")
        }
}

process IntersectVCF {
    publishDir "${params.output_dir}/variants/contaminants", pattern: "*contaminants.vcf", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(sample_id), val(file_id), path(variants), path(synthetics)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}_contaminants.vcf")

    script:
        """
        bedtools intersect -header -a ${variants} -b ${synthetics}  > ${file_id}_contaminants.vcf
        """
}

process AnnotateVCF {
    publishDir "${params.output_dir}/variants/annotated", pattern: "*annotated.vcf", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(sample_id), val(file_id), path(vcf)

    output:
        tuple val(sample_id), val(file_id), path("${file_id}_annotated.vcf")

    script:
        """
        annotate_vcf.py ${vcf} ${file_id}_annotated.vcf
        """
}

process FreebayesContaminants {
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_huge_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), val(file_id), path(input_bam_file), path(input_bai_file), path(reference), path(reference_idx)
        tuple val(roi_sample_id), val(roi_file_id), path(roi)

    output:
        tuple val(sample_id), path(file_id), path("${input_bam_file.SimpleName}.multiallelic.vcf")

    script:
        """
        freebayes \
        -f $reference \
        --min-alternate-fraction 0 \
        --min-alternate-count 1 \
        --min-base-quality 0 \
        --vcf ${input_bam_file.SimpleName}.multiallelic.vcf \
        $input_bam_file
        """
}

process SeparateMultiallelicVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(vcf)

    output:
        tuple val(X), path("${vcf.SimpleName}.vcf")

    script:
        """
        bcftools norm -m -any $vcf > ${vcf.SimpleName}.vcf
        """
}

process FindRegionOfInterest{
    //  publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'
    
    input:
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)

    output:
        tuple val(sample_id), val(file_id), path("${sample_id}_roi.bed")

    script:
        """
        samtools depth $bam_in | awk '\$3>${params.roi_detection.min_depth}' | awk '{print \$1"\t"\$2"\t"\$2 + 1}' | bedtools merge -d ${params.roi_detection.max_distance} -i /dev/stdin > ${sample_id}_roi.bed
        """
}


// REPORTING
process CountFastqInfo{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
        val(mode)

    output:
        path ("${sample_id}_${mode}_read_count.txt")
        path ("${sample_id}_${mode}_base_count.txt")
        path ("${sample_id}_${mode}_overview.txt")

    
    script:
        """
        seqkit stats -T $fastq | tee ${sample_id}_${mode}_overview.txt | awk 'BEGIN{fs = "\t"} { sum+=\$4} END{print sum}' > ${sample_id}_${mode}_read_count.txt
        cat ${sample_id}_${mode}_overview.txt | tr -d \\, | awk 'BEGIN{fs = "\t"} { sum+=\$5} END{print sum}' > ${sample_id}_${mode}_base_count.txt
        """
}

process PlotFastqsQUalAndLength {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq)
    val grep_pattern
    val mode
    val tab_name

    output:
    tuple val(sample_id), path("${sample_id}_${mode}_histograms.json")

    script:
    """
    plot_fastq_histograms.py ${grep_pattern} ${sample_id}_${mode}_histograms.html ${tab_name} ${params.priority_limit} \
    2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${sample_id}_${mode}_histograms.json
    """
}

process PlotReadStructure {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_cpu_medium'

    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sample_id), val(file_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${bam.simpleName}_read_structure.json")

    script:
    // This takes a lot of RAM when the sequencing summary is big!
    """
    samtools sort -n -o tmp_readname_sorted_${bam.simpleName}.bam ${bam}
    plot_read_structure_donut.py tmp_readname_sorted_${bam.simpleName}.bam \
    ${bam.simpleName}_read_structure.json \
    ${params.priority_limit} \
    2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${bam.simpleName}_read_structure.json
    """
}

process PlotVcf {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sample_id), val(file_id), path(vcf)

    output:
    tuple val(sample_id), path("${vcf.simpleName}.json")

    script:
    // This takes a lot of RAM when the sequencing summary is big!
    """
    plot_vcf.py ${vcf} ${vcf.simpleName}.html ${params.priority_limit} \
    2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${vcf.simpleName}.json
    """
}

process PasteVariantTable {
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(file_id), path(vcf_file)

    output:
    tuple val(sample_id), path("${vcf_file.simpleName}_table.json")

    script:
    """
    write_variants_table.py ${vcf_file} ${vcf_file.simpleName}_table.json 'Variant table' ${params.priority_limit} \
    2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${vcf_file.simpleName}_table.json
    """
}

process PasteContaminantTable {
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id), path("${sample_id}_contaminant_table.json")

    script:
    """
    write_contaminants_table.py ${vcf_file} ${sample_id}_contaminant_table.json 'Contaminants' ${params.priority_limit} \
    2> >(tee -a ${sample_id}_error.txt >&2) || catch_plotting_errors.sh ${sample_id}_error.txt ${sample_id}_contaminant_table.json
    """
}

process PlotQScores {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(split_pileup), path(consensus_pileup)

    output:
    tuple val(sample_id), path("${consensus_pileup.simpleName}.json")

    script:
    """
    plot_bam_accuracy.py ${split_pileup} ${consensus_pileup} ${consensus_pileup.simpleName}.html ${params.priority_limit} \
    2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${consensus_pileup.simpleName}.json
    """
}

process PlotMetadataStats {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_huge_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sample_id), path(jsons)

    output:
    tuple val(sample_id), path("${sample_id}_metadata_plots.json")

    script:
    """
    plot_metadata.py . ${sample_id}_metadata_plots.html ${params.priority_limit} --subsample_size ${params.metadata.subsample_size} \
    2> >(tee -a ${sample_id}_error.txt >&2) || catch_plotting_errors.sh ${sample_id}_error.txt "${sample_id}_metadata_plots.json"
    """
}

process PlotReport {
    publishDir "${params.output_dir}", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
    tuple val(sample_id), path(jsons)

    output:
    tuple val(sample_id), path("${sample_id}_report.html")

    script:
    """
    generate_report.py ${sample_id} '${params}' ${workflow.manifest.version} ${params.priority_limit}
    """
}

process SamtoolsQuickcheck{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)

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
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)
    
    output:
        tuple val(sample_id), path("${bam_in.SimpleName}.flagstats_metadata.json")

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
        tuple val(sample_id), val(file_id), path(bam_in), path(bai_in)
    
    output:
        tuple val(sample_id), path("${bam_in.SimpleName}.idxstats_metadata.json")

    script:
        """
        TOTALREFREADS=\$(samtools idxstats $bam_in | grep -v "^BB" | grep -v "^*" | awk -F' ' '{sum+=\$3;} END{print sum;}')
        echo \$TOTALREFREADS
        TOTALREFREADS=\$TOTALREFREADS jq -n '{additional_info:{"total_reference_mapping_reads":env.TOTALREFREADS,}}' > ${bam_in.SimpleName}.idxstats_metadata.json

        """
}

process CountNonBackboneVariants{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), path(vcf)

    output:
        tuple val(sample_id), path("${vcf.SimpleName}.variants.metadata.json")

    script:
        """
        VARS=\$(cat $vcf | grep -v "#" | wc -l)
        NONBBVARS=\$(cat $vcf | grep -v "#" | grep -v "^BB" | wc -l)
        VARS=\$VARS NONBB=\$NONBBVARS jq -n '{additional_info:{"variants_found":env.VARS, "variants_found_non_backbone": env.NONBB}}' > ${vcf.SimpleName}.variants.metadata.json
        """
}