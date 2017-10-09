#!/usr/local/bin/Rscript --verbose
library(optparse)
source('FastIndep.R')
source('kingraph.R')

cat(as.character(Sys.time()),"\n\n", file = stderr())

# Setup command line options  ##################################################

option_list <- list(
  make_option(c("-t","--phi"), type = "numeric", default = "0.0442", 
              help = "phi threshold to infer independence [default= %default]",
              metavar = "option"),
  make_option(c("-n","--nreq"), type = "numeric", default = "2", 
              help = "number of ivs requested to fastindep [default= %default]",
              metavar = "option"),
  make_option(c("-m","--nread"), type = "numeric", default = "2", 
              help = "number of ivs to read from fastindep [default= %default]",
              metavar = "option"),
  make_option(c("-v", "--vertexinfo"), type = "character", 
              help="tab delimited file, minimum 2 columns: name collection",
              metavar="option"),
  make_option(c("-r", "--resultsdir"), type = "character", default = "fastivs",
              help="results directory [default= %default]",
              metavar="option")
); 

opt_parser <- OptionParser(usage = "%prog [options] file",
                           option_list=option_list);

# Initialization   #############################################################
args <- parse_args2(opt_parser)
opts <- args$options

input <- args$args

if (length(input) == 0) {
  stop("No input file given.\n",print_help(opt_parser))
} else{
  print(args$args)
  file <- args$args[1] # Normalize path?
}


if (!is.null(opts$vertexinfo)) {
  vertex.info <- read.table(file=opts$vertexinfo, header=TRUE)
  if (!all(c('name','collection') %in% colnames(vertex.info))){
    stop(paste('Either `name` or `collection` headers missing from',
               opts$vertexinfo))
  }
} else{
  vertex.info <-data.frame(name = character(0),
                           collection = character(0)) 
}

results.dir <- file.path(normalizePath(dirname(opts$resultsdir)),
                         basename(opts$resultsdir))

dir.create(results.dir)

thresh <- opts$phi
n <- opts$nreq
m <- opts$nread

file.base <- basename(tools::file_path_sans_ext(file))
output.base <- file.path(results.dir,
                         paste(file.base,thresh, sep="."))

# Find  maximum inpendent vertex set ###########################################
# (with fast indep greedy by default)

phi <- read.phi(input) 
phi.mat <- paste(output.base,"matrix.tab",sep = ".")

write.table(adj.from.phi(phi)$matrix,
            file = phi.mat,
            quote = FALSE, sep = "\t")

fastindep.out <- paste(output.base,"fastindep",sep = ".")
fastindep.log <- paste(output.base,"fastindep.log",sep = ".")

ivs.run <- run.fastindep(thresh = thresh,
                   n = n,
                   input = phi.mat,
                   output = fastindep.out,
                   log = fastindep.log)

ivs.out <- read.fastindep(fastindep.out, max = m)

# Minimum Maximal independent set ##############################################
# hints to parental lines
adj <- adj.from.phi(phi)
g <- kinship.graph(adj, remove.unrelated = FALSE)
is.ivs <- as.factor(V(g)$name %in% ivs.out$best)
independent <- data.frame(name = V(g)$name,
                          ivs = is.ivs)
independent$color.border <- c(NA,'black')[is.ivs]

if (nrow(vertex.info) > 0){
  vertex.info <- merge(vertex.info,independent, by='name')
  as.type = 'collection'
} else{
  vertex.info <- independent
  as.type = NULL
}

# Add vertex set  info to graph ################################################
g <- add.vertex.info(g,vertex.info)

# Text Output ##################################################################
# print max.ivs.list
# print summary table
# file.base phi.threshold	dup.fr n	n.independent	n.redundant

write.table(ivs.out$best,
            paste(output.base, ".fastivs.list", sep=""),
            row.names = FALSE,
            col.names = FALSE,
            quote =FALSE, sep="\t")

write.table(V(g)$name[degree(g) == 0],
            paste(output.base, ".singlets.list", sep=""),
            row.names = FALSE,
            col.names = FALSE,
            quote =FALSE, sep="\t")

# Plot Output ##################################################################
# plot adjacency matrix # print to pdf?
# plot html mst

# adj.heatmap(adj, pdf.base = paste(output.base, "matrix", sep="."))

legend.arg <- list(addNodes = data.frame(label = c("HapMap","CIAT",
                                                   "independent set"), 
                                         shape = c("diamond","dot","dot"),
                                         color.border = c(NA,NA,"black")))

g.vis <- vis.kingraph(color.components(mst(g)), legend.arg = legend.arg)

visSave(g.vis[[2]],
        paste(output.base,"mst","html", sep="."),
        selfcontained = TRUE, background = "lightgrey")

system2("rm", args = c("-r", paste(results.dir,"*_files", sep = "/")))
