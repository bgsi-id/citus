nextflow.enable.dsl = 2

include { CITUS } from './workflows/citus'

workflow {
    CITUS()
}
