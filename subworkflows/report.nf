
nextflow.enable.dsl=2

// Run the report generator on the jsons provided

include {
    GenerateHtmlReport
} from '../modules/reporting'

workflow  Report{
    take:
        json_reads
        json_globals
        vcf
    main:
        GenerateHtmlReport(json_reads,
            json_globals,
            vcf
        )
    emit:
        GenerateHtmlReport.output
}
