// perform transcript assembly & quant with bambu across multiple samples as a merged analysis
include { BAMBU_ASSEMBLY as BAMBU_MERGE       } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_MERGE_QUANT } from '../../../modules/local/bambu/assembly/main'

workflow MULTISAMPLE_TRANSCRIPT_QUANT {
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
    BAMBU_MERGE(
        merge_ch,
        yieldsize,
        NDR,
        fasta,
        gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_MERGE.out.versions)
    // run bambu in quantification mode
    if (quant) {
        multi_gtf_ch = BAMBU_MERGE.out.transcriptome.map { meta, dir ->
            def merge_gtf = file("${dir}/extended_annotations.gtf")
            merge_gtf
        }
        multi_gtf_ch.view()
        BAMBU_MERGE_QUANT(
            bam_ch,
            yieldsize,
            NDR,
            fasta,
            multi_gtf_ch,
        )
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
