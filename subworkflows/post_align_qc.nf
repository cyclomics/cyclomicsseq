#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Run some simple checks when alignment is done

include {
    SamtoolsQuickcheck
    SamtoolsFlagstats
} from "./modules/samtools"


workflow  PostAlignmentQC{
    take:
        final_alignment_bam  
    main:
        SamtoolsQuickcheck(final_alignment_bam)
        SamtoolsQuickcheck.out.view()
        SamtoolsFlagstats(final_alignment_bam)
        SamtoolsFlagstats.out.view()

}
