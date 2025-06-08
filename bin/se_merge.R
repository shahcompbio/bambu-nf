#!/usr/bin/env Rscript
library(bambu)
library(optparse)
## load in commandline args
option_list = list(
  make_option("--se", type="character", default=NULL, 
              help="summarized experiments", metavar="character"),
  make_option("--out_dir", type="character", default=NULL, 
              help="output directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# create vector of summarized experiment objects
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}
rdata_paths <- strsplit(opt$se, ",")[[1]]
### combine objects
se <- loadRData(rdata_paths[1])
for (i in 2:length(rdata_paths)){
  se1 <- loadRData(rdata_paths[i])
  se <- cbind(se, se1)
}
# write bambu outputs
writeBambuOutput(se, path=opt$out_dir)
se_out <- sprintf("%s/merged_se.RData", opt$out_dir)
save(se, file = se_out)