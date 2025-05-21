#!/usr/bin/env Rscript
library(bambu)
## load in commandline args
args = commandArgs(trailingOnly=TRUE)
sample.bam <- strsplit(grep('--bam*', args, value = TRUE), split = '=')[[1]][[2]]
yieldSize <- strsplit(grep('--yieldsize*', args, value = TRUE), split = '=')[[1]][[2]]
##prep annotations ...
fa.file <- strsplit(grep('--ref_genome*', args, value = TRUE), split = '=')[[1]][[2]]
gtf.file <- strsplit(grep('--ref_gtf*', args, value = TRUE), split = '=')[[1]][[2]]
# This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
annotations <- prepareAnnotations(gtf.file)
## create RC files
se.discoveryOnly <- bambu(reads = sample.bam, annotations = annotations, ncore = 1, genome = fa.file,
            yieldSize = yieldSize, lowMemory=TRUE, fusionMode=FALSE,
            verbose=TRUE, discovery=TRUE, quant=FALSE, rcOutDir = "./")
