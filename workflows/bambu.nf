/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_READS                                      } from '../subworkflows/local/preprocess_reads/main'
include { SINGLE_TRANSCRIPT_QUANT                               } from '../subworkflows/local/single_transcript_quant/main'
include { MULTISAMPLE_TRANSCRIPT_QUANT                          } from '../subworkflows/local/multisample_transcript_quant/main'
include { MULTISAMPLE_TRANSCRIPT_QUANT as MULTISAMPLE_FIXED_NDR } from '../subworkflows/local/multisample_transcript_quant/main'
include { MULTIQC                                               } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                } from '../subworkflows/local/utils_nfcore_bambu-nf_pipeline'
// modules for merge workflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAMBU_NF {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run samtools view to filter bam files for reads aligned to accessory chromosomes
    //
    if (!params.skip_preprocessing) {
        input_ch = ch_samplesheet.map { meta, bam, bai, rds -> tuple(meta, bam, bai) }
        PREPROCESS_READS(input_ch, params.filter_reads, params.filter_acc_reads)
        rc_ch = PREPROCESS_READS.out.reads
        bam_ch = PREPROCESS_READS.out.bam
        ch_versions = ch_versions.mix(PREPROCESS_READS.out.versions)
    }
    else {
        rc_ch = ch_samplesheet.map { meta, bam, bai, rds -> tuple(meta, rds) }
        bam_ch = ch_samplesheet.map { meta, bam, bai, rds -> tuple(meta, bam) }
    }
    // perform assembly & quantification with bambu
    if (params.single_sample) {
        SINGLE_TRANSCRIPT_QUANT(rc_ch, params.recommended_NDR, params.yieldsize, params.fasta, params.gtf, params.NDR)
        ch_versions = ch_versions.mix(SINGLE_TRANSCRIPT_QUANT.out.versions)
    }
    if (!params.skip_multisample) {
        merge_ch = rc_ch
            .collect { meta, rds -> rds }
            .map { rds -> [["id": "merge"], rds] }
        // run at recommended NDR
        if (params.recommended_NDR) {
            MULTISAMPLE_TRANSCRIPT_QUANT(
                merge_ch,
                bam_ch,
                [],
                params.yieldsize,
                params.fasta,
                params.gtf,
                params.quantification,
            )
            ch_versions = ch_versions.mix(MULTISAMPLE_TRANSCRIPT_QUANT.out.versions)
        }
        if (params.NDR != null) {
            MULTISAMPLE_FIXED_NDR(
                merge_ch,
                bam_ch,
                params.NDR,
                params.yieldsize,
                params.fasta,
                params.gtf,
                params.quantification,
            )
            ch_versions = ch_versions.mix(MULTISAMPLE_FIXED_NDR.out.versions)
        }
    }
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'bambu-nf_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
