// QC samples
include { NANOPLOT } from '../../../modules/nf-core/nanoplot/main'

workflow QC {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    NANOPLOT(ch_bam)
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] })
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)

    emit:
    nanostats = NANOPLOT.out.txt // channel: [ val(meta), [ txt ] ]
    versions  = ch_versions // channel: [ versions.yml ]
    multiqc   = ch_multiqc_files // channel: [ path(multiqc) ]
}
