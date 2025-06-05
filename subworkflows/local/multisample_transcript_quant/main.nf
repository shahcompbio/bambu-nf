// perform transcript assembly & quant with bambu across multiple samples as a merged analysis
include { BAMBU_ASSEMBLY as BAMBU_MERGE     } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_MERGE_NDR } from '../../../modules/local/bambu/assembly/main'
include { BAMBU_ASSEMBLY as BAMBU_MERGE_QUANT     } from '../../../modules/local/bambu/assembly/main'

workflow MULTISAMPLE_TRANSCRIPT_QUANT {
    take:
    merge_ch        // channel: [ val(meta), [ rds ] ]
    bam_ch       // channel: [ val(meta), path(bam) ]
    recommended_NDR // boolean; use recommended NDR for bambu assembly
    yieldsize       // integer; number of reads to process with bambu
    fasta           // ref genome
    gtf             // ref gtf
    NDR             // float; fixed NDR value for bambu assembly

    main:

    ch_versions = Channel.empty()
    // run bambu merge at different NDRs
    if (recommended_NDR) {
        BAMBU_MERGE(
            merge_ch,
            yieldsize,
            [],
            fasta,
            gtf,
        )
        ch_versions = ch_versions.mix(BAMBU_MERGE.out.versions)
    }
    // run at fixed NDR
    BAMBU_MERGE_NDR(
        merge_ch,
        yieldsize,
        NDR,
        fasta,
        gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_MERGE_NDR.out.versions)
    // run bambu in quantification mode
    multi_gtf_ch = BAMBU_MERGE
                    .out.transcriptome
                    .map { meta, dir -> 
                    def gtf = file("${dir}/extended_annotations.gtf")
                    gtf}
    multi_gtf_ch.view()
    BAMBU_MERGE_QUANT(
        bam_ch,
        yieldsize,
        [],
        fasta,
        multi_gtf_ch,
    )
    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
