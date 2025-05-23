#!/usr/bin/env Rscript
library(bambu)
library(optparse)
## load in commandline args
option_list = list(
  make_option("--rds", type="character", default=NULL, 
              help="bambu read classes", metavar="character"),
  make_option("--yieldsize", type="character", default=1e5, 
              help="limits number of reads processed at once", metavar="character"),
  make_option("--NDR", type="numeric", default=NULL, 
              help="bambu NDR; modulates FDR", metavar="character"),
  make_option("--ref_genome", type="character", default=NULL, 
              help="reference genome", metavar="character"),
  make_option("--ref_gtf", type="character", default=NULL, 
              help="reference gtf", metavar="character"),
  make_option("--out_dir", type="character", default=NULL, 
              help="output dir", metavar="character")
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
writeBambuOutput(se, path = opt$out_dir)
#
## let's also save the summarized experiment object to play with later ...
se_out <- sprintf("%s/se.RData", opt$out_dir)
save(se, file = se_out)