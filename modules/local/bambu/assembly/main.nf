// transcript assembly with bambu
process BAMBU_ASSEMBLY {
    tag "${meta.id}"
    // for testing purposes
    label 'process_single'
    publishDir "${params.outdir}/${meta.id}/transcriptome_NDR_${task.ext.args}", mode: 'copy', overwrite: true

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/bambu:3.10.0"

    input:
    tuple val(meta), path(rds)
    val yieldsize
    path ref_genome
    path ref_gtf

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    transcript_assembly.R \\ 
        --rds=${rds} \\ 
        --yieldsize=${yieldsize} \\ 
        --ref_genome=${ref_genome} \\ 
        --ref_gtf=${ref_gtf} \\ 
        --NDR=${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bambu: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
