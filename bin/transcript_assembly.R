#!/usr/bin/env Rscript
library(bambu)
## load in commandline args
args = commandArgs(trailingOnly=TRUE)
rds <- strsplit(grep('--rds*', args, value = TRUE), split = '=')[[1]][[2]]
yieldSize <- strsplit(grep('--yieldsize*', args, value = TRUE), split = '=')[[1]][[2]]
##prep annotations ...
fa.file <- strsplit(grep('--ref_genome*', args, value = TRUE), split = '=')[[1]][[2]]
gtf.file <- strsplit(grep('--ref_gtf*', args, value = TRUE), split = '=')[[1]][[2]]
# This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
annotations <- prepareAnnotations(gtf.file)
####################
### run with recommended NDR
#############
## make output for multisample using recommended NDR
multi_out <- sprintf("%s/transcriptome_NDR_default/", bambu_out)
create_directory(multi_out)
ncores <- length(rds)
##trying to set baselineFDR vs NDR here to be more permissive ...
## baselineFDR defaults to 0.1 when NDR is not selected ...
## try with model training to begin with as Andre Sim recommended ...
se.NDR_default <- bambu(reads = rds, annotations = annotations, ncore = ncores,
            genome = fa.file, verbose=TRUE, yieldSize = yieldSize, lowMemory=TRUE,
            fusionMode=FALSE)

## write outputs to gtf and expression level files
writeBambuOutput(se.NDR_default, path = multi_out)
#
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.NDR_default.RData", multi_out)
save(se.NDR_default, file = se_out)
### now let's test different NDR thresholds
## start by doing transcript discovery
newAnnotations <- bambu(reads = rds, annotations = annotations, ncore = ncores,
                        genome = fa.file, NDR = 1, verbose=TRUE,
                        yieldSize = yieldSize, lowMemory=TRUE, quant = FALSE)
##################
### do quantification on all new transcripts (NDR = 1) ....
###########
se.NDR_1 <- bambu(reads = rds, annotations = newAnnotations, ncore = ncores,
                  genome = fa.file, verbose=TRUE, yieldSize = yieldSize,
                  lowMemory=TRUE, NDR = 1, discovery = FALSE)
multi_out <- sprintf("%s/transcriptome_NDR_1/", bambu_out)
create_directory(multi_out)
## write outputs to gtf and expression level files
writeBambuOutput(se.NDR_1, path = multi_out)
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.NDR_1.RData", multi_out)
save(se.NDR_1, file = se_out)
####################
## do quantification on NDR < 0.1 transcripts ...
##################
## make a filtered annotation
annotations.filtered <- newAnnotations[(!is.na(mcols(newAnnotations)$NDR) & mcols(newAnnotations)$NDR<0.1) | is.na(mcols(newAnnotations)$NDR)]
se.NDR_0.1 <- bambu(reads = rds, annotations = annotations.filtered, ncore = ncores,
                    genome = fa.file, verbose=TRUE, yieldSize = yieldSize,
                    lowMemory=TRUE, NDR = 1, discovery = FALSE)
multi_out <- sprintf("%s/transcriptome_NDR_0.1/", bambu_out)
create_directory(multi_out)
## write outputs to gtf and expression level files
writeBambuOutput(se.NDR_0.1, path = multi_out)
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.NDR_0.1.RData", multi_out)
save(se.NDR_0.1, file = se_out)