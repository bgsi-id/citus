nextflow.enable.dsl = 2

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'
// def summary_params = paramsSummaryMap(workflow)

include { DUMPVERSION }                             from '../modules/local/dumpversion.nf'
include { PB_GERMLINE }                             from '../modules/local/parabricks.nf'

include { FASTQC }                                  from '../modules/nf-core/fastqc/main.nf'
include { SAMTOOLS_CONVERT }                        from '../modules/nf-core/samtools/convert/main.nf'
include { SAMTOOLS_INDEX }                          from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_STATS }                          from '../modules/nf-core/samtools/stats/main.nf'
include { MOSDEPTH }                                from '../modules/nf-core/mosdepth/main.nf'
include { TABIX_BGZIP }                             from '../modules/nf-core/tabix/bgzip/main'   
include { TABIX_TABIX }                             from '../modules/nf-core/tabix/tabix/main'   
include { BCFTOOLS_STATS }                          from '../modules/nf-core/bcftools/stats/main.nf'
include { VCFTOOLS as VCFTOOLS_TSTV_COUNT }         from '../modules/nf-core/vcftools/main.nf'
include { VCFTOOLS as VCFTOOLS_TSTV_QUAL }          from '../modules/nf-core/vcftools/main.nf'
include { VCFTOOLS as VCFTOOLS_SUMMARY }            from '../modules/nf-core/vcftools/main.nf'
include { MULTIQC }                                 from '../modules/nf-core/multiqc/main.nf'

workflow CITUS {
  if (params.samplesheet_input != 'NO_FILE') {
    fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [ [id: it['ID']], [ it['R1'], it['R2'] ] ] }
  } else {
    fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  multiqc_config = "${workflow.projectDir}/assets/multiqc_config.yml"
  description = "${workflow.projectDir}/assets/methods_description_template.yml"

  mqc_config  = Channel.fromPath( multiqc_config, checkIfExists: true )
  ch_multiqc_custom_methods_description = file( description, checkIfExists: true  )

  fasta       = Channel.fromPath( params.fasta ).first()
  fai         = Channel.fromPath( "${params.fasta}.fai" ).first()
  target      = Channel.fromPath( params.target ).first()
  bwa         = Channel.fromPath( params.bwa_index ).collect()
  
  reports     = Channel.empty()
  versions    = Channel.empty()

  main:
    FASTQC(fastq)
    reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
    versions = versions.mix(FASTQC.out.versions.first())
    
    PB_GERMLINE(fastq, fasta, bwa, params.known_site)
    reports = reports.mix(PB_GERMLINE.out.recal.collect{ meta, report -> report })
    versions = versions.mix(PB_GERMLINE.out.versions)
    
    SAMTOOLS_INDEX(
      PB_GERMLINE.out.bam
    )

    bam_bai = PB_GERMLINE.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)

    SAMTOOLS_CONVERT(
      bam_bai,
      fasta.map{ it -> [ [ id:'fasta' ], it ] },
      fai.map{ it -> [ [ id:'fai' ], it ] }
    )

    SAMTOOLS_STATS(
      bam_bai,
      fasta.map{ it -> [ [ id:'fasta' ], it ] }
    )
    reports = reports.mix(SAMTOOLS_STATS.out.stats.collect{ meta, report -> report })
    versions = versions.mix(SAMTOOLS_STATS.out.versions)

    MOSDEPTH(bam_bai.combine(target), fasta.map{ it -> [ [ id:'fasta' ], it ] })
    reports = reports.mix(MOSDEPTH.out.global_txt.collect{ meta, report -> report })
    reports = reports.mix(MOSDEPTH.out.regions_txt.collect{ meta, report -> report })
    versions = versions.mix(MOSDEPTH.out.versions)

    TABIX_BGZIP(
      PB_GERMLINE.out.vcf
    )

    TABIX_TABIX(
      TABIX_BGZIP.out.output
    )
    versions = versions.mix(TABIX_TABIX.out.versions)

    vcf_tbi = TABIX_BGZIP.out.output.join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    BCFTOOLS_STATS(
      vcf_tbi,
      target.map{ it -> [ [ id:'target' ], it ] },
      [[:],[]],
      [[:],[]],
      [[:],[]],
      [[:],[]],
    )
    VCFTOOLS_TSTV_COUNT(
      PB_GERMLINE.out.vcf,
      target,
      []
    )
    VCFTOOLS_TSTV_QUAL(
      PB_GERMLINE.out.vcf,
      target,
      []
    )
    VCFTOOLS_SUMMARY(
      PB_GERMLINE.out.vcf,
      target,
      []
    )
    reports = reports.mix(BCFTOOLS_STATS.out.stats.collect{ meta, stats -> stats })
    reports = reports.mix(VCFTOOLS_TSTV_COUNT.out.tstv_count.collect{ meta, stats -> stats })
    reports = reports.mix(VCFTOOLS_TSTV_QUAL.out.tstv_qual.collect{ meta, stats -> stats })
    reports = reports.mix(VCFTOOLS_SUMMARY.out.filter_summary.collect{ meta, stats -> stats })

    versions = versions.mix(BCFTOOLS_STATS.out.versions)
    versions = versions.mix(VCFTOOLS_SUMMARY.out.versions)

    version_yaml = Channel.empty()
    DUMPVERSION(versions.unique().collectFile(name: 'collated_versions.yml'))
    version_yaml = DUMPVERSION.out.mqc_yml.collect()

    // workflow_summary    = WorkflowCitus.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCitus.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    multiqc_files = Channel.empty()
    multiqc_files = multiqc_files.mix(version_yaml)
    // multiqc_files = multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    multiqc_files = multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    multiqc_files = multiqc_files.mix(reports.collect().ifEmpty([]))
    MULTIQC(multiqc_files.collect(), mqc_config.collect().ifEmpty([]), [], [])
}