
nextflow.enable.dsl=2

// Run the report generator on the jsons provided

include {
    JqAddDepthToJson
} from "./modules/jq"

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
        full_json = JqAddDepthToJson(json_globals, depth_json)

        if (params.control_vcf){
            control = Channel.fromPath(params.control_vcf, checkIfExists: true)
        report = GenerateHtmlReportWithControl(json_reads,
            full_json,
            vcf,
            control
            )
        }
        else {
            report = GenerateHtmlReport(json_reads,
            full_json,
            vcf
            )
        }
    emit:
        report
}
