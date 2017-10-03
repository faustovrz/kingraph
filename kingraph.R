library(reshape2)
library(igraph)
library(Matrix)
library(visNetwork)
library(htmltools)
library(RColorBrewer)
library(gplots)


phi.palette <- c("blue","purple","green","yellow","orange","red")
phi.breaks <- c(0,0.0442,0.0884,0.177,0.354,0.45,0.5)
phi.legend <- c("Clones >= 0.45 ",
                "1st deg inbred [0.354,0.45)",
                "1st deg [0.1777,0.354)",
                "2nd deg [0.0884,0.177)", 
                "3rd deg [0.0442,0.0884)",
                "Unrelated < 0.0442",
                "No het")

read.phi <- function(file.name="relatedness2", header=TRUE){
  # ASSUME standard vcftools relatedness2 output table
  # INDV1	INDV2	N_AaAa	N_AAaa	N1_Aa	N2_Aa	RELATEDNESS_PHI
  skip <- 0
  if(header){ skip <- 1} 
  df <-read.table(file=file.name,header=header, na.strings = c("-nan","-inf"))
  colnames(df) <- c("INDV1","INDV2","N_HetHet","N_AAaa","N1_Aa","N2_Aa","PHI")
  if (nrow(df[is.na(df$PHI),]) > 0) {
    print(df[is.na(df$PHI),])
    warning("WARNING: Non numeric kinship coefficients found and removed!")
  }
  
  df <- df[!is.na(df$PHI),]
  df
  #phi <- df[as.numeric(df$INDV1) < as.numeric(df$INDV2),c("INDV1","INDV2","PHI")]
}


adj.from.phi <- function(phi, thresh = -100, negtozero = TRUE){
  edges <- phi[phi$PHI >= thresh,c("INDV1","INDV2","PHI")]
  edges$weight <- 10 * (1 - edges$PHI)
  # weight similarity more than distance
  # replace negative values with 0
  if (negtozero == TRUE & any(edges$PHI < 0)  ){
  edges[which(edges$PHI < 0), ]$weight <- 0
  }
  df <- dcast(edges, INDV1~INDV2, value.var="PHI")
  m <- data.matrix(df[,-1])
  if (negtozero == TRUE ) { m[m < 0] <- 0 } # redundant with above conditional
  m[is.na(m)] <- 0
  rownames(m) <- colnames(m)
  list(edges = edges, matrix = m)
}

write.kinship.matrix <- function(adj, file="phi_matrix.tab"){
  write.table(adj$matrix,file="phi_square", quote=FALSE, sep="\t")
}


adj.heatmap <- function(adj, pdf.base="phi.heatmap.pdf",phi.thresh = 0.0442,symm=TRUE){
  if(is.matrix(adj)){matrix <- adj}
  else if(is.matrix(adj$matrix)){ matrix <- adj$matrix}
  dims <- ceiling(dim(matrix)/100)*10
  pdf(paste(pdf.base,".", phi.thresh,".pdf",sep=""), height=dims[1], width=dims[2])
  heatmap.2(matrix, col=phi.palette, breaks = phi.breaks, key=FALSE, trace="none",
            symm=symm,symkey=F,symbreaks=F)
  
  legend(xpd=TRUE, 0,1,bty = "n",
         legend=phi.legend,
         title = "Kinship Coefficient",
         cex=1, fill=c(rev(phi.palette),"white"))
  dev.off()
}



vis.kingraph <- function(g,legend.arg=NULL, bipartite = FALSE, remove.phi=NULL){
  
  if (!is.null(remove.phi)) { 
    g <- remove.edges.by.phi(g,phi = remove.phi)
  }
  
  v <- visIgraph(g)  %>% visOptions(  width = "1280px", height="960px", nodesIdSelection = TRUE, highlightNearest = TRUE)
  if(!is.list(legend.arg) & ( is.character(V(g)$collection) |  is.factor(V(g)$collection) )) {
    legend.arg <- list(addNodes = data.frame(label = c("HapMap","CIAT","independent set"), 
                                             shape = c("diamond","dot","dot"),
                                             color.border = c(NA,NA,"black")))
  }
  ledges <- data.frame(color = rev(phi.palette)[1:5], 
                       label = phi.legend[1:5], font.align = "bottom",
                       width = 5,
                       arrows = "")
  legend.arg$addEdges=ledges
  legend.arg$useGroups=FALSE
  v <- do.call(visLegend, c(list(v),legend.arg))
  if (bipartite == TRUE){
   v <- v %>% visIgraphLayout("layout_as_bipartite") 
  }
  
  browsable(
    tagList(
      tags$head(
        tags$style('div.vis-network{background-color: lightgrey;}')  
      ),
      v
    )
  )
}


kinship.graph <- function(adj, thresh = 0.0442, deg = NULL, remove.unrelated =TRUE){
  df <- adj$edges[as.numeric(adj$edges$INDV2) > as.numeric(adj$edges$INDV1),]
  
  if (remove.unrelated) { df <- df[df$PHI >= thresh,] }
  
  if (any(df$PHI >= thresh)) {
  edges <- as.matrix( df[c("INDV1", "INDV2")])
  g <- graph_from_edgelist(edges, directed = FALSE)
  
  V(g)$title <- V(g)$name
  E(g)$weight <- df$weight
  E(g)$PHI <- df$PHI
  
  # Remove edges corresponding to unrelated pairs phi < 0.0442 if any
  # TODO: Add argument to remove singleton nodes
  # from the kinship network and modify thresh accordingly
  # by default I will be adding them
  if (any(E(g)$PHI < 0.0442)){
  g <- delete.edges(g, E(g)[E(g)$PHI < 0.0442])
  }
  g <- delete.edges(g, E(g)[E(g)$PHI < thresh])
  
  E(g)$color <- phi.palette[cut(E(g)$PHI, phi.breaks)]
  # Sort by name
  # TODO: maybe get this its own function!!!!
  g <- permute(g, invPerm((order(V(g)$name))))
  V(g)$title <- V(g)$name
  
  g
  } else {
    make_empty_graph(directed = FALSE)
  }
}


phi.cc.ivs <- function(cc,phi, scale = 50){
if (ecount(cc) ==0) { 
    return( list(ivs.size=1,
                 independent=V(cc),
                 # simplify this
                 redundant= V(cc)[!(V(cc)$name %in% independent$name)], 
                 g=cc) )
}
  
if (!(ecount(cc) * vcount(cc) <= scale**2)) {
  stop(paste("Highly dense graph, try with a sparser one.\n", 
             round(ecount(cc)/vcount(cc), digits=1), " edges per vertex.",
             vcount(cc)), " vertices")
  } else {
  print(paste(round(ecount(cc)/vcount(cc), digits=1), " edges per vertex.", vcount(cc), " vertices"))

  # I am picking the first solution arbitrarily.
  # Good enough!
  # TODO: filter by solutions containning a specific set of vertices 
  independent <- ivs(cc, min=ivs_size(cc))[[1]]  
  list(ivs.size=ivs_size(cc),
       independent=independent,
       redundant=V(cc)[!(V(cc)$name %in% independent$name)],
       g=cc)
  }
  
}

phi.ivs <- function(g,phi,scale = 50){
  #comp_vcount <- unlist(lapply(decompose(g),vcount))
  ivs.list <- lapply(decompose(g),phi.cc.ivs,phi,scale=scale)
  independent.name  <-names( unlist (sapply(ivs.list, "[[", "independent")))
  redundant.name <-names( unlist (sapply(ivs.list, "[[", "redundant")))
  ivs.size <- sum(unlist (sapply(ivs.list, "[[", "ivs.size")))
  V(g)$color[V(g)$name %in% independent.name] <- "red"
  
  list(ivs.size=ivs.size, 
       independent=V(g)[V(g)$name %in% independent.name], 
       redundant=V(g)[!(V(g)$name %in% independent.name)],
       g=g,
       phi = phi[!(phi$INDV1 %in% redundant.name) & !(phi$INDV2 %in% redundant.name),])
}


kin.wash <- function(ivs, thresh=0.45){
  # adj.from.phi(ivs$phi)$edges
  g <- kinship.graph(adj.from.phi(ivs$phi), thresh = thresh)
  phi.ivs(kinship.graph(adj.from.phi(ivs$phi), thresh = thresh),
                   ivs$phi)
  # g <- kinship.graph(adj.from.phi(ivs1$phi), thresh = 0.0442)
  # adj.heatmap(adj.from.phi(ivs1$phi), pdf.base=paste("phi",thresh))
}

kin.stepIVS <- function(phi, upper=0.45, lower=0.0442, nv = 10){
  phi.ud <- phi[as.numeric(phi$INDV2) > as.numeric(phi$INDV1),]
  phi.related <- phi.ud$PHI[phi.ud$PHI>=lower]
  if (max(phi.related) < upper) { 
      upper <- max(phi.related)
    }
  step <- (nv^2)/(4*(length(phi.related)))
  q <- quantile( phi.related, seq(0,1,by=step))
  q <- rev( c(q[q > lower & q < upper]))
  print(paste("upper :", upper, "q1 :", q[1]))
  q[length(q)] <- min(c(lower, q[length(q)] ))

  # with the quantile steps a solvable IVS problem is set
  # quartz()
  # plot(ecdf(phi.related),xlim=c(0,0.5))
  kin.g <-  kinship.graph(adj.from.phi(phi), remove.unrelated = FALSE)
  step.g <- kinship.graph(adj.from.phi(phi), thresh = q[1])
  step.ivs <- phi.ivs(step.g, phi)
  step.redundant <- step.ivs$redundant
  print(paste("phi: ",q[1],". ","Initial ivs step out of", length(q)))

  i <- 0
  for (phi_step in q[-1]) {
    if (any(step.ivs$phi$PHI >= phi_step)) { # this validation makes no sense in 
                                             # the phi matrix the diagonal is always 0.5
    step.ivs <- kin.wash(step.ivs, thresh=phi_step)
    step.redundant <- c( step.redundant, step.ivs$redundant)
    i <- i+1
    print(paste("phi: ",phi_step,". ", "ivs step ", i , "out of",length(q) - 1) )
    } else{ 
    i <- i+1
    print(paste("phi: ",phi_step,". ", "ivs step ", i , "out of", length(q) - 1))
    next } 
  }
 step.ivs$phi <- phi[ !(phi$INDV1 %in% names(step.redundant)) & !(phi$INDV2 %in% names(step.redundant)),]
 step.ivs$independent <- V(kin.g)[!(V(kin.g)$name %in% names(step.redundant))]
 step.ivs$redundant <- V(kin.g)[V(kin.g)$name %in% names(step.redundant)]
 step.ivs$g <- color.components(kin.g)
 V(step.ivs$g)$color.border[V(step.ivs$g)$name %in% names(step.ivs$independent)] <- "black"
 step.ivs
}

color.components<- function(g){
c <- components(g)
mst <- mst(g)
V(mst)[match(names(c$membership), V(g)$name)]$group  <-  c$membership
V(g)[match(names(c$membership), V(g)$name)]$group  <-  c$membership
g
}


add.collection <- function(g,col.df){
  # add collection (CIAT vs HapMap) data.
  # TODO: generalize to add any dataframe to the graph
  # for exmple cultivated vs wild, geographic region: AFR, CAR, MES, SAM
  g.col <- merge(x=data.frame(name=V(g)$name), y=col.df, by="name")
  shp <- c("dot","diamond")
  type <- c(TRUE,FALSE)
  V(g)$collection <- shp[g.col$collection]
  V(g)$type <-  type[g.col$collection]
  V(g)$shape <- shp[g.col$collection]
  V(g)$title <-  V(g)$name
  g
}

# Duplicate Network
find.clones <- function(phi,col= data.frame()){
  g <- kinship.graph(adj.from.phi(phi), thresh = 0.45)
  if(nrow(col.df)>0){
    g <-add.collection(g,col.df)
  }
  clone.ivs <- phi.ivs(g, phi)
  
  V(g)$color.border[V(clone.ivs$g)$name %in% names(clone.ivs$independent)] <- "black"

  if(length(col.df) == 0){
  legend.arg <- list(addNodes = data.frame(label = c("independent set"), 
                                             shape = c("dot"),
                                             color.border = c("black")))
  }else{ 
  g <- add.collection(g,col.df)
  legend.arg <- list(addNodes = data.frame(label = c("CIAT","HapMap","independent set"), 
                                             shape = c("dot","diamond","dot"),
                                             color.border = c(NA,NA,"black")))
    
  }  
  list(g = g, 
       redundant = clone.ivs$redundant,
       independent  = clone.ivs$independent,
       phi = phi[ !(phi$INDV1 %in% names(clone.ivs$redundant)) & !(phi$INDV2 %in%  names(clone.ivs$redundant)),],
       legend.arg = legend.arg)
}


kinship.analysis <- function(phi.file,col.df=data.frame(),ivs.thresh=0.177, results.dir="results"){
  phi <- read.phi(phi.file) 
  file.base <-tools::file_path_sans_ext(phi.file)
  file.base <-basename(file.base)
  
  # maybe getwd() is redundant
  
  dir.create(file.path(getwd(),results.dir))
  
  output.base <- file.path(getwd(),results.dir, paste(file.base,ivs.thresh, sep="."))
  adj <- adj.from.phi(phi)

  # kinship network
  
  g <- kinship.graph(adj,remove.unrelated=FALSE, thresh = -200)
  g.thresh <-  kinship.graph(adj,remove.unrelated=FALSE, thresh = ivs.thresh)
  # percentage of related individuals
  vertex.count <-  vcount(g)
  unrelated.count <- sum(degree(g.thresh)==0)
  related.pct <- 1 - unrelated.count/vertex.count
  connectivity <- ecount(g.thresh)/(vertex.count*(vertex.count-1)/2)
  adj.heatmap(adj, pdf.base = paste(output.base,".matrix", sep=""),phi.thresh = -100)
  
  dup<- find.clones(phi,col.df)
  phi045 <- phi$PHI >0.45
  clone.independent <- phi$INDV1  %in% names(dup$redundant)
  clone.redundant <- phi$INDV2  %in% names(dup$independent)
  
  write.table(phi[clone.independent & clone.redundant & phi045,],
              paste(output.base,".dup.phi", sep=""),
              row.names = FALSE,quote=FALSE, sep="\t")
  nr.df <- data.frame(nr = unique(c(as.character(dup$phi$INDV1), 
                                    as.character(dup$phi$INDV2))))
  write.table(nr.df,
              paste(output.base,".nodup.list", sep=""),
              row.names = FALSE, col.names= FALSE, quote=FALSE, sep="\t")
  
  
  visSave(vis.kingraph(color.components(dup$g))[[2]],
          paste(output.base,".nodup.html", sep=""),
          selfcontained = TRUE, background = "lightgrey")
  if(length(col.df) > 0){
    if(length(levels(as.factor(V(dup$g)$type)))==2 ){
      col1 <- which(col.df$collection==levels(col.df$collection)[1])
      col2 <- which(col.df$collection==levels(col.df$collection)[2])
      adj.heatmap(adj$matrix[col2,col1], symm=FALSE,
                  pdf.base = paste(output.base,".bp.matrix", sep=""),
                  phi.thresh = -100)
      
      bp <- remove.degree.0(kin.bipartite.graph(dup$g))
      visSave(vis.kingraph(bp,bipartite=TRUE)[[2]],
              paste(output.base,".dup.bp.html", sep=""),
            selfcontained = TRUE, background = "lightgrey")
    }
  }
  adj.heatmap(adj.from.phi(dup$phi),
              pdf.base = paste(output.base,".nodup.matrix", sep=""),
              phi.thresh = -100)
  
  # g <- add.collection(g,col)
  # vis.kingraph(mst(g), legend.arg = legend.arg)
  
  
  #kinship connectivity
  
  # Independent Vertex Set
  # kin.ivs<- kin.stepIVS(phi, upper=0.45, lower=ivs.thresh, nv = 10)
  # beginning at dup$phi removes duplicates from analysis hereafter
  kin.ivs<- kin.stepIVS(dup$phi, upper=0.45, lower=ivs.thresh, nv = 10)
  kin.ivs$g <- delete_edges(kin.ivs$g,which(E(kin.ivs$g)$PHI < ivs.thresh))
  ivs.df <- data.frame(id =names(kin.ivs$independent))
  
  write.table(ivs.df,
              paste(output.base,".ivs.list", sep=""),
              row.names = FALSE, col.names= FALSE, quote=FALSE, sep="\t")
  
  if(length(col.df) > 0){
    kin.ivs$g <- add.collection(kin.ivs$g, col.df)
    if(length(levels(as.factor(V(kin.ivs$g)$type))) == 2 ){
      legend.arg <- list(addNodes = data.frame(label = c("CIAT","HapMap","independent set"), 
                                               shape = c("dot","diamond","dot"),
                                               color.border = c(NA,NA,"black")))
    } else{
      legend.arg <- list(addNodes = data.frame(label = c("independent set"), 
                                               shape = c("dot"),
                                               color.border = c("black")))
      }
    
    }else{
      legend.arg <- list(addNodes = data.frame(label = c("independent set"), 
                                               shape = c("dot"),
                                               color.border = c("black")))
    }
    
  kin.graph <- vis.kingraph(color.components(mst(kin.ivs$g)),legend.arg = legend.arg)[[2]]

  visSave(kin.graph, paste(output.base,".mst.html", sep=""),
          selfcontained = TRUE, background = "lightgrey")
  
  if(length(col.df) > 0)
    if(length(levels(as.factor(V(kin.ivs$g)$type))) ==2){
      bp <- kin.bipartite.graph(kin.ivs$g)
      visSave(vis.kingraph(bp,bipartite = TRUE )[[2]], paste(output.base,".bp.html", sep=""),
                           selfcontained = TRUE, background = "lightgrey")
      
     visSave(vis.kingraph(mst(bp))[[2]], paste(output.base,".mst.bp.html", sep=""),
              selfcontained = TRUE, background = "lightgrey")
      
      visSave(vis.kingraph(mst(bp),bipartite = TRUE)[[2]], paste(output.base,".mst.bp.suriyama.html", sep=""),
                           selfcontained = TRUE, background = "lightgrey")
    
    }
  
  
  kin.adj <- adj.from.phi(kin.ivs$phi)
  
  
  adj.heatmap(adj.from.phi(kin.ivs$phi),
              pdf.base = paste(output.base,".ivs.matrix", sep=""),
              phi.thresh = -100)
  
  data.frame ( file.base=file.base,
               phi.thresh=ivs.thresh,
               related.fr = related.pct,
               connectivity =connectivity,
               with.dup = vcount(dup$g),
               dup.independent = length(dup$independent),
               dup.redundant = length(dup$redundant),
               dup.fr = vcount(dup$g)/vcount(kin.ivs$g),
               n = vcount(kin.ivs$g),
               n.independent = length(kin.ivs$independent),
               n.redundant = length(kin.ivs$redundant))
}

kin.bipartite.graph <-function(g, type.attr = "type", remove.unrelated=FALSE){
  #check bipartite
  type <- as.factor(get.vertex.attribute(g, type.attr))
  type <- as.factor(get.vertex.attribute(g, "collection")) 
  if (length(levels(type)) != 2){
    stop(paste(type.attr, " is not a bipartite attribute: ", levels(type)))
     
  }
    type <- c(TRUE,FALSE)[type]
    g <- set_vertex_attr(g, "type", V(g), type)

  part1 <- V(g)$name[which(V(g)$type==FALSE)]
  part2 <- V(g)$name[which(V(g)$type==TRUE)]
  bp.combination <- expand.grid(part1,part2)
  colnames(bp.combination) <-c("V1","V2")
  
  bp.edges <-c(t(merge(as.data.frame(get.edgelist(g)),
                       bp.combination,
                       by=c("V1","V2"))))
  
  bp.edges.ids <-get.edge.ids(g, bp.edges)
  
  bp <- delete.edges(g, E(g)[!(E(g) %in% E(g)[bp.edges.ids])])
  if(remove.unrelated == TRUE){
    bp <- delete.vertices(bp,which(degree(bp)<1))
  }
  bp
}

remove.edges.by.phi<- function(g, phi=0.0442){
  g <- delete.edges(g, E(g)[E(g)$PHI < phi])
}

remove.degree.0<- function(g){
  g <- delete.vertices(g, which(degree(g)==0))
}

