nextflow.enable.dsl = 2


include { validateParameters; paramsHelp } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def String command = "nextflow run"
log.info logo + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)

// Print help message if needed
if (params.help) {
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log)

include { CITUS } from './workflows/citus'

workflow {
    CITUS()
}
