// perform transcript assembly & quant with bambu across multiple samples as a merged analysis
include { BAMBU_ASSEMBLY as BAMBU       } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_QUANT } from '../../../modules/local/bambu/assembly/main'

workflow TRANSCRIPT_QUANT {
    take:
    merge_ch  // channel: [ val(meta), [ rds ] ]
    bam_ch    // channel: [ val(meta), path(bam) ]
    NDR       // float; empty value will use recommended NDR
    yieldsize // integer; number of reads to process with bambu
    fasta     // ref genome
    gtf       // ref gtf
    quant     // boolean; true if running in quantification mode

    main:

    ch_versions = Channel.empty()
    // run bambu merge
    BAMBU(
        merge_ch,
        yieldsize,
        NDR,
        fasta,
        gtf,
    )
    ch_versions = ch_versions.mix(BAMBU.out.versions)
    // run bambu in quantification mode
    if (quant) {
        multi_gtf_ch = BAMBU.out.transcriptome.map { meta, dir ->
            def merge_gtf = file("${dir}/extended_annotations.gtf")
            merge_gtf
        }
        multi_gtf_ch.view()
        BAMBU_QUANT(
            bam_ch,
            yieldsize,
            NDR,
            fasta,
            multi_gtf_ch,
        )
    }
    // emit summarized experiment for further processing
    if (quant) {
        ch_se = BAMBU_QUANT.out.se
    }
    else {
        ch_se = BAMBU.out.se
    }
    ch_se.view()
    // add NDR to metamap
    ch_se = ch_se.map { meta, se -> [meta + ["NDR": NDR], se] }

    emit:
    se       = ch_se // channel: [ val(meta), path(se) ]
    versions = ch_versions // channel: [ versions.yml ]
}
