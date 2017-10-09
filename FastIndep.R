library(stringr)

# FastIndep  wrapper ###########################################################
# Example
# indep <- run.fastindep(thresh = 0.0442,
#                       n = 100,
#                       input = 'phi/WILDLAC.HapMap.745.esculenta.phi.matrix.txt',
#                       output = 'phi/WILDLAC.HapMap.745.esculenta.indep',
#                      log = 'phi/rfastindep.log')

run.fastindep <- function(exe = 'fastindep', thresh=0.0442,
                          n=2, input="", output="", log=""){
  # TODO: add tmp file ?
  
  args <- paste(c('-t', '-n', '-i', '-o'),
                c(thresh, n, input, output))
  
  return(system2(exe, args = args, stdout = log))
}

# FastIndep Parser ###########################################################

# Example
# fast.ivs <- read.fastindep("phi/WILDLAC.HapMap.745.esculenta.indep", 100)
# hist(fast.ivs$set.size)
# fast.ivs$greedy

read.fastindep <-function(file,max=2){
  #stop if max >2
  if (max >2){
    stop(paste("read at least 2 sets from",file))
  }
  # TODO: Add validation of file size, 1GB?
  f <- file(description=file, open="r")
  stop <- FALSE
  set.list <- list()
  set.size <- c()
  n.set <- 0
  result.log <- c()

while(!stop) {
    line <- readLines(f, n = 1)
    if(length(line) == 0) {
      stop <- TRUE
      break()
    }
    
    line <- gsub("^\\s","",line, perl = TRUE)

    if(line == ""){
      next
    } else if (grepl("Set Size =",line)){
      n.set <- n.set + 1
      if (n.set <= max){
        line <- gsub("Set Size =\\s+\\d+\\s+","",line, perl = TRUE)
        set <- strsplit(line,"\\s+", perl =TRUE)
        set.list[n.set] <- set
      }
    } else if(grepl("Map Size",line)){
      size.count <- str_extract_all(line, "[0-9]+", simplify = TRUE)
      size.count <- as.numeric(size.count)
      set.size<-c(set.size,rep(size.count[1],size.count[2]))
    } else {
      result.log <- c(result.log,line)
    }
    
  }
  close(f)
  greedy.size <- length(set.list[[1]])
  random.size <- length(set.list[[2]])
  max.size <- max(greedy.size,random.size)
  
  is.max <-c(greedy.size,random.size) == max.size
  list(n.set= n.set, set.list= set.list,set.size= set.size,
       best = set.list[[which(is.max)]],
       greedy = set.list[[1]],
       random = set.list[[2]],
       log = result.log)
}


