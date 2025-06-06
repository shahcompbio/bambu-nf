// merge summarized experiments
process SEMERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/bambu:3.10.0beta"

    input:
    tuple val(meta), path(se, arity: '1..*')

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.RData"), emit: se
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def se_list = se.join(',')
    """
#!/usr/bin/env Rscript
## load RData function
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}
rdata_paths <- strsplit(${se_list}, ",")[[1]]
### combine objects
se <- loadRData(rdata_paths[1])
for (i in 2:length(rdata_paths)){
  se1 <- loadRData(rdata_paths[i])
  se <- cbind(se, se1)
}
### save experiment
save(se, file = "all_se.RData")
## write out a small YAML with the R version:
ver <- R.Version()\$version.string
out <- sprintf('"%s":\\n    R: "%s"\\n', "${task.process}", ver)
cat(out, file="versions.yml")
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
        semerge: \$(semerge --version)
    END_VERSIONS
    """
}
