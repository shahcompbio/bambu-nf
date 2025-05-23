#!/usr/bin/env Rscript
library(bambu)
library(optparse)
## load in commandline args
option_list = list(
  make_option(c("--rds"), type="character", default=NULL, 
              help="bambu read classes", metavar="character"),
  make_option(c("--yieldsize"), type="character", default="1000000", 
              help="limits number of reads processed at once", metavar="character"),
  make_option(c("--NDR"), type="character", default=NULL, 
              help="bambu NDR; modulates FDR", metavar="character"),
  make_option(c("--ref_genome"), type="character", default=NULL, 
              help="reference genome", metavar="character"),
  make_option(c("--ref_gtf"), type="character", default=NULL, 
              help="reference gtf", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
annotations <- prepareAnnotations(opt$ref_gtf)
####################
se <- bambu(reads = opt$rds, annotations = annotations, ncore = 1,
            genome = opt$ref_genome, NDR = opt$NDR, verbose=TRUE, 
            yieldSize = opt$yieldsize, lowMemory=TRUE)
## write outputs to gtf and expression level files
writeBambuOutput(se, path ="./")
#
## let's also save the summarized experiment object to play with later ...
save(se, file = "./se.RData")