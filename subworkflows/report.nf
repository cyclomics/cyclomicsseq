
nextflow.enable.dsl=2

// Run the report generator on the jsons provided

include {
    GenerateHtmlReport
    GenerateHtmlReportWithControl
} from './modules/reporting'

workflow  Report{
    take:
        json_reads
        json_globals
        vcf
    main:
        // main has a if statement as workaround for conditional input parameters
        if (params.control_vcf){
            control = Channel.fromPath(params.control_vcf, checkIfExists: true)
        report = GenerateHtmlReportWithControl(json_reads,
            json_globals,
            vcf,
            control
            )
        }
        else {
            report = GenerateHtmlReport(json_reads,
            json_globals,
            vcf
            )
        }
    emit:
        report
}
