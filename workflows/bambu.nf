/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_VIEW               } from '../modules/nf-core/samtools/view/main'
include { BAMBU_READCLASSES           } from '../modules/local/bambu/readclasses/main'
include { BAMBU_ASSEMBLY ; BAMBU_ASSEMBLY as BAMBU_NDR } from '../modules/local/bambu/assembly/main'
include { BAMBU_FILTER                } from '../modules/local/bambu/filter/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_bambu-nf_pipeline'

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
    if (params.filter_acc_reads) {
        ch_filtered = SAMTOOLS_VIEW(
            ch_samplesheet,
            [[], []],
            [],
            "bai",
        )
        ch_bam = ch_filtered.bam
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
    }
    else {
        ch_bam = ch_samplesheet.map { meta, bam, bai -> tuple(meta, bam) }
    }
    // create read classes with bambu
    rc_ch = BAMBU_READCLASSES(
        ch_bam,
        params.yieldsize,
        params.fasta,
        params.gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_READCLASSES.out.versions)
    // perform assembly & quantification with bambu
    if (params.recommended_NDR) {
        ch_bambu_default = BAMBU_ASSEMBLY(
            rc_ch.rds,
            params.yieldsize,
            [],
            params.fasta,
            params.gtf,
        )
        ch_versions = ch_versions.mix(BAMBU_ASSEMBLY.out.versions)
    }
    // run at fixed NDR
    ch_bambu_ndr = BAMBU_NDR(
        rc_ch.rds,
        params.yieldsize,
        params.NDR,
        params.fasta,
        params.gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_NDR.out.versions)
    // filter for detected transcripts
    BAMBU_FILTER(ch_bambu_ndr.se)
    ch_versions = ch_versions.mix(BAMBU_FILTER.out.versions)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
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
