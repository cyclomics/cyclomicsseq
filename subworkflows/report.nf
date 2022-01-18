
nextflow.enable.dsl=2

// Run the report generator on the jsons provided

include {
    GenerateHtmlReport
} from '../modules/reporting'

workflow  Report{
    take:
        jsons
    emit:
        GenerateHtmlReport.output
    main:
        GenerateHtmlReport(jsons)
}
