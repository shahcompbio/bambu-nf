/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_READS                           } from '../subworkflows/local/preprocess_reads/main'
include { BAMBU_ASSEMBLY as BAMBU } from '../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_MERGE } from '../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_MERGE_QUANT } from '../modules/local/bambu/assembly/main'
include { TRANSCRIPT_QUANT ; TRANSCRIPT_QUANT as TRANSCRIPT_QUANT_MERGE } from '../subworkflows/local/transcript_quant/main'
include { SEMERGE                                    } from '../modules/local/semerge/main'
include { BAMBU_FILTER                               } from '../modules/local/bambu/filter/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_bambu-nf_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAMBU_NF {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // shorten param names called by workflows for readability
    def yield = params.yieldsize
    def genome = params.fasta
    def gtf = params.gtf
    def merge_quant = params.multisample_quant
    // begin workflow
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
    // make an NDR channel
    if (params.recommended_NDR && params.NDR != null) {
        ch_NDR = channel.of([], params.NDR)
    }
    else if (params.recommended_NDR) {
        ch_NDR = channel.of([])
    }
    else {
        ch_NDR = channel.of(params.NDR)
    }
    // closure to handle empty NDR
        def padNDR = { tup ->
        tup.size() < 3
            ? [ tup[0], tup[1], [] ]
            : tup
    }
    ref_gtf_ch = channel.of(params.gtf)

    // run single sample assembly & quant with read classes
    if (params.single_sample) {
        // combine read classes and NDR channels
        bambu_input_ch = rc_ch
                            .combine(ch_NDR)
                            .map(padNDR)
                            .combine(ref_gtf_ch)
        BAMBU(bambu_input_ch, yield, genome)
        ch_versions = ch_versions.mix(BAMBU.out.versions)
    }
    // run multisample mode; assembly then quant
    if (!params.skip_multisample) {
        merge_ch = rc_ch
            .collect { meta, rds -> rds }
            .map { rds -> [["id": "merge"], rds] }
        merge_input_ch = merge_ch
                            .combine(ch_NDR)
                            .map(padNDR)
                            .combine(ref_gtf_ch)
        BAMBU_MERGE(merge_input_ch, yield, genome)
        ch_versions = ch_versions.mix(BAMBU_MERGE.out.versions)
        // run multisample quantification
        if (merge_quant) {
            // run bambu in quantification mode with merged gtf
            ext_gtf_ch = BAMBU_MERGE.out.gtf.map{meta, NDR, ext_gtf -> tuple(NDR, ext_gtf)}.view()
            merge_quant_ch = bam_ch.combine(ext_gtf_ch)
            BAMBU_MERGE_QUANT(merge_quant_ch,yield,genome)
            // collect all summarized experiments for each NDR
            //BAMBU_MERGE_QUANT.out.se.view()
            // se_ch = TRANSCRIPT_QUANT_MERGE.out.se
            //     .map { meta, se -> tuple(meta.NDR, se) }
            //     .groupTuple()
            // se_ch.view()
        }
    }
    // filter for detected transcripts
    // all_se_ch = merge_se_ch.mix(TRANSCRIPT_QUANT.out.se)
    // BAMBU_FILTER(all_se_ch)
    // ch_versions = ch_versions.mix(BAMBU_FILTER.out.versions)
    // collect versions
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
