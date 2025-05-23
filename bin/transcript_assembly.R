#!/usr/bin/env Rscript
library(bambu)
## load in commandline args
args = commandArgs(trailingOnly=TRUE)
rds <- strsplit(grep('--rds*', args, value = TRUE), split = '=')[[1]][[2]]
yieldSize <- strsplit(grep('--yieldsize*', args, value = TRUE), split = '=')[[1]][[2]]
NDR <- strsplit(grep('--NDR*', args, value = TRUE), split = '=')[[1]][[2]]
fa.file <- strsplit(grep('--ref_genome*', args, value = TRUE), split = '=')[[1]][[2]]
gtf.file <- strsplit(grep('--ref_gtf*', args, value = TRUE), split = '=')[[1]][[2]]
# This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
annotations <- prepareAnnotations(gtf.file)
####################
se <- bambu(reads = rds, annotations = annotations, ncore = 1,
            genome = fa.file, NDR = NDR, verbose=TRUE, 
            yieldSize = yieldSize, lowMemory=TRUE)
## write outputs to gtf and expression level files
writeBambuOutput(se, path ="./")
#
## let's also save the summarized experiment object to play with later ...
save(se, file = "./se.RData")