// perform transcript assembly & quant with bambu across multiple samples as a merged analysis
include { BAMBU_ASSEMBLY as BAMBU_MERGE     } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_MERGE_NDR } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_QUANT     } from '../../../modules/local/bambu/assembly/main'

workflow MULTISAMPLE_TRANSCRIPT_QUANT {
    take:
    merge_ch        // channel: [ val(meta), [ rds ] ]
    recommended_NDR // boolean; use recommended NDR for bambu assembly
    yieldsize       // integer; number of reads to process with bambu
    fasta           // ref genome
    gtf             // ref gtf
    NDR             // float; fixed NDR value for bambu assembly

    main:

    ch_versions = Channel.empty()

    // merge transcriptomes across multiple samples
    merge_ch.view()
    // run bambu merge at different NDRs
    if (recommended_NDR) {
        ch_bambu_merge = BAMBU_MERGE(
            merge_ch,
            yieldsize,
            [],
            fasta,
            gtf,
        )
        ch_versions = ch_versions.mix(BAMBU_MERGE.out.versions)
    }
    // run at fixed NDR
    ch_bambu_ndr = BAMBU_MERGE_NDR(
        merge_ch,
        yieldsize,
        NDR,
        fasta,
        gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_MERGE_NDR.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
