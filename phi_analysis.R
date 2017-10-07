#!/usr/local/bin/Rscript --verbose
library(optparse)
source('kingraph.R')

cat(as.character(Sys.time()),"\n\n", file = stderr())

phi.help <- "comma separated string defining phi values: [default= %default]"

option_list <- list(
  make_option("--phi", type = "character", default = "0.177,0.0442", 
              help= phi.help,
              metavar = "character"),
  make_option(c("-c", "--col"), type = "character", 
              help="collection table with 2 columns: name collection",
              metavar="character")
); 

opt_parser <- OptionParser(usage = "%prog [options] file1 file2 ...",
                           option_list=option_list);

args <- parse_args2(opt_parser)
opts <- args$options

input <- args$args

if (length(input) == 0) {
  stop("No input file given.\n",print_help(opt_parser))
}

kin.coeff <- as.numeric(unlist(strsplit(opts$phi, ",")))

if (!is.null(opts$col)) {
  col.df <- read.table(file=opts$col, header=TRUE)
  if (any(!(c('name','collection') %in% colnames(col.df)))){
    stop(paste('Either `name` or `collection` headers missing from',
             opts$col))
  }
} else{
  col.df <-data.frame(name = character(0), collection = character(0)) 
}


dir.create("results")

columns <- c("file.base", "phi.thresh", "related.fr", "connectivity",
             "with.dup", "dup.independent", "dup.redundant", "dup.fr",
             "n", "n.independent", "n.redundant")

write.table(t(columns), file= "results/summary.tab",row.names = FALSE,
            quote=FALSE, sep="\t", col.names =FALSE)

for (file in input){
  for(coeff in kin.coeff){
    kin.summary <- kinship.analysis(file, 
                                    col= col.df,
                                    ivs.thresh = coeff,
                                    results.dir= "results")
    write.table(kin.summary, file= "results/summary.tab",
                quote=FALSE, sep="\t", append = TRUE,
                row.names=FALSE, col.names = FALSE)
  }
}

system("rm -r results/*_files")

