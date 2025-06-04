// perform transcript assembly & quant with bambu across individual samples
// no merged analysis

include { BAMBU_ASSEMBLY ; BAMBU_ASSEMBLY as BAMBU_NDR } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_FILTER                } from '../../../modules/local/bambu/filter/main'
workflow SINGLE_TRANSCRIPT_QUANT {
    take:
    rc_ch           // channel: [ val(meta), [ rds ] ]
    recommended_NDR // boolean; use recommended NDR for bambu assembly
    yieldsize       // integer; number of reads to process with bambu
    fasta           // ref genome
    gtf             // ref gtf
    NDR             // float; fixed NDR value for bambu assembly

    main:

    ch_versions = Channel.empty()
    // perform assembly & quantification with bambu
    if (recommended_NDR) {
        ch_bambu_default = BAMBU_ASSEMBLY(
            rc_ch,
            yieldsize,
            [],
            fasta,
            gtf,
        )
        ch_versions = ch_versions.mix(BAMBU_ASSEMBLY.out.versions)
    }
    // run at fixed NDR
    ch_bambu_ndr = BAMBU_NDR(
        rc_ch,
        yieldsize,
        NDR,
        fasta,
        gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_NDR.out.versions)
    // filter for detected transcripts
    BAMBU_FILTER(ch_bambu_ndr.se)
    ch_versions = ch_versions.mix(BAMBU_FILTER.out.versions)

    emit:
    versions         = ch_versions // channel: [ versions.yml ]
}

