// preprocess reads via filtering and read class creation with bambu

include { SAMTOOLS_VIEW ; SAMTOOLS_VIEW as FILTER_READS } from '../../../modules/nf-core/samtools/view/main'
include { BAMBU_READCLASSES             } from '../../../modules/local/bambu/readclasses/main'

workflow PREPROCESS_READS {
    take:
    ch_samplesheet   // channel: [ val(meta), [ bam, bai, rcFile ] ]
    filter_reads     // boolean; filter reads on mapq and read length
    filter_acc_reads // boolean; filter reads on accessory chromosomes

    main:

    ch_versions = Channel.empty()

    if (filter_reads) {
        ch_filtered = FILTER_READS(
            ch_samplesheet,
            [[], []],
            [],
            "bai",
        )
        ch_bam = ch_filtered.bam
        ch_versions = ch_versions.mix(FILTER_READS.out.versions.first())
    }
    else if (filter_acc_reads) {
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
    BAMBU_READCLASSES(
        ch_bam,
        params.yieldsize,
        params.fasta,
        params.gtf,
    )
    ch_versions = ch_versions.mix(BAMBU_READCLASSES.out.versions)

    emit:
    bam      = ch_bam // channel: [ val(meta), path(bam) ]
    reads    = BAMBU_READCLASSES.out.rds // channel: [ val(meta), [ rcFile ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
