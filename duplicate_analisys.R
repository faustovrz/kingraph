#!/usr/local/bin/Rscript --verbose
library(optparse)
source('kingraph.R')

cat(as.character(Sys.time()),"\n\n", file = stderr())

# Setup command line options  ##################################################

option_list <- list(
  make_option(c("-t","--phi"), type = "numeric", default = "0.45", 
              help = "phi threshold to infer duplication [default= %default]",
              metavar = "option"),
  make_option(c("-v", "--vertexinfo"), type = "character", 
              help="tab delimited file, minimum 2 columns: name collection",
              metavar="option"),
  make_option(c("-r", "--resultsdir"), type = "character", default = "duplicate",
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

# Analysis per file  ###########################################################

for (file in input){
  # Setup IO ###################################################################
  file.base <- basename(tools::file_path_sans_ext(file))
  output.base <- file.path(results.dir,
                           paste(file.base,thresh, sep="."))
  
  # Kinship network ############################################################
  
  phi <- read.phi(file) 
  adj <- adj.from.phi(phi)
  g <-kinship.graph(adj, remove.unrelated = FALSE)
  dup<- find.dup(phi,vertex.info = vertex.info)

  dup.vis <- vis.kingraph(color.components(dup$g),
                          legend.arg = dup$legend.arg)
  
  is.dup <- phi$PHI > thresh
  INDV1.is.redundant <- phi$INDV1  %in% names(dup$redundant)
  INDV2.is.independent <- phi$INDV2  %in% names(dup$independent)
  
  nodup <- sort(unique(c(as.character(dup$phi$INDV1), 
                         as.character(dup$phi$INDV2))))
  
  if(is.bipartite(dup$g)){
    comp <- compare.bipart(dup$g,adj)
    dup.bp.vis <- vis.kingraph(comp$g, bipartite = TRUE,
                              legend.arg = dup$legend.arg)
  }

  stats <- kin.stats(g)
  dup.stats <- kin.stats(dup$g)
  
  dup.summary  <- list(file.base = file.base,
                       phi.thresh = thresh,
                       n = vcount(g),
                       related.fr =  stats$related.fr,
                       connectivity = stats$connectivity,
                       with.dup =  vcount(dup$g),
                       dup.independent = length(dup$independent),
                       dup.redundant = length(dup$redundant),
                       dup.fr =  vcount(dup$g)/vcount(g))
  
  # Text output ################################################################
  write.table(phi[ is.dup & INDV1.is.redundant & INDV2.is.independent, ],
              paste(output.base,".dup.phi", sep=""),
              row.names = FALSE, quote=FALSE, sep="\t")
  
  write.table(dup$phi,
              paste(output.base,".nodup.phi", sep=""),
              row.names = FALSE, quote=FALSE, sep="\t")
  cat(nodup,
      file =  paste(output.base,"nodup","list", sep ="."),
      sep = "\n")
  
  write.table(dup.summary,
              paste(output.base,".dup.summary", sep=""),
              row.names = FALSE, quote=FALSE, sep="\t")
  
  # Plot output ################################################################
  
  adj.heatmap(adj,pdf.base = paste(output.base,"matrix", sep="."))
  
  adj.heatmap(adj.from.phi(dup$phi),
              pdf.base = paste(output.base, "nodup.matrix", sep="."))
  
  visSave(dup.vis[[2]],
          paste(output.base,"dup","html", sep="."),
          selfcontained = TRUE, background = "lightgrey")
  
  if(is.bipartite(dup$g)){
  adj.heatmap(t(comp$adj), symm=FALSE,
              pdf.base = paste(output.base,"bp.matrix", sep="."))
    
  visSave(dup.bp.vis[[2]],
          paste(output.base,"dup.bp.html", sep="."),
          selfcontained = TRUE, background = "lightgrey")
  }
}

system2("rm", args = c("-r", paste(results.dir,"*_files", sep = "/")))


