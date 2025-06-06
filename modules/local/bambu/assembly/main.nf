// transcript assembly with bambu
// NOTE: be warying of errant spaces in the script section! optparse will be unhappy
process BAMBU_ASSEMBLY {
    tag "${meta.id}_NDR_${NDR ?: 'DEFAULT'}"
    cpus { rds.size() > 20 ? 20 : rds.size() }
    // for testing purposes
    label 'process_high_memory'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', overwrite: true

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/bambu:3.10.0beta"

    input:
    tuple val(meta), path(rds, arity: '1..*'), val(NDR), path(ref_gtf)
    val yieldsize
    path ref_genome

    output:
    tuple val(meta), val(NDR), path("*/se.RData"), emit: se
    tuple val(meta), val(NDR), path("transcriptome_*"), emit: transcriptome
    tuple val(meta), val(NDR), path("*/extended_annotations.gtf"), emit: gtf, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def NDR_args = NDR ? "--NDR=${NDR}" : ""
    def out_dir = NDR ? "transcriptome_NDR_${NDR}" : "transcriptome_NDR_DEFAULT"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def rds_list = rds.join(',')
    """
    mkdir -p ${out_dir}
    transcript_assembly.R \\
        --rds=${rds_list} \\
        --yieldsize=${yieldsize} \\
        --ref_genome=${ref_genome} \\
        --ref_gtf=${ref_gtf} \\
        --out_dir=${out_dir} \\
        --ncore=${task.cpus} \\
        ${NDR_args} ${args}
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
