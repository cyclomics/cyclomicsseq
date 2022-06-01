
nextflow.enable.dsl=2

// Run the report generator on the jsons provided

include {
    AddDepthToJson
} from "./modules/bin"

include {
    GenerateHtmlReport
    GenerateHtmlReportWithControl
} from './modules/reporting'


workflow  Report{
    take:
        json_reads
        json_globals
        vcf
        depth_json
    main:
        // main has a if statement as workaround for conditional input parameters
        // full_json = AddDepthToJson(json_globals, depth_json)
        full_json = json_globals
        
        if (params.control_vcf){
            println('Generating report with control vcf')
            control = Channel.fromPath(params.control_vcf, checkIfExists: true)
        report = GenerateHtmlReportWithControl(json_reads,
            full_json,
            vcf,
            control
            )
        }
        else {
            println('Generating report without control vcf')
            report = GenerateHtmlReport(json_reads,
            full_json,
            vcf
            )
        }
    emit:
        report
}
