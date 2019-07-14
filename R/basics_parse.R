#' Parse haps file
#'
#' This function uses fread to parse a .haps file, which is part of the haps/sample output file format of SHAPEIT2.
#'
#' @param filename string. Filename
#' @param ... Any other parameters for fread
#' @return Returns a data table.
#' @examples
#' read.haps("example.haps")
#'

read.haps <- function(filename, ...){
  haps <- fread(filename, ...)
  return(haps)
}

#' Parse sample file
#'
#' This function uses fread to parse a .sample file, which is part of the haps/sample output file format of SHAPEIT2.
#'
#' @param filename string. Filename
#' @param ... Any other parameters for fread
#' @return Returns a data table.
#' @examples
#' read.haps("example.sample")
#'
read.sample <- function(filename, ...){
  sam <- fread(filename, ...)
  return(sam)
}

#' Parse mut file
#'
#' This function uses fread to parse a .mut file, which is part of the anc/mut output file format of SHAPEIT2.
#'
#' @param filename string. Filename
#' @param CHR int. Chromosome index.
#' @param ... Any other parameters for fread
#' @return Returns a data table.
#' @examples
#' read.haps("example.mut")
#'
read.mut <- function(filename, CHR = NA, ...){
  mut <- fread(filename, ...)
  mut <- mut[,1:(ncol(mut)-1)]
  if(!is.na(CHR)) mut   <- cbind(CHR = CHR, mut)

  if(ncol(mut) < 11){
    warning("Specified mut file has too few columns.")
  }

  return(mut)
}

#' Parse coal file
#'
#' This function parses a .coal file, which stores coalescence rates inferred by Relate
#'
#' @param filename string. Filename
#' @return Returns a data frame.
#' @examples
#' read.coal("example.coal")
#'
read.coal <- function(filename, ...){

  groups <- as.matrix(read.table(filename, nrow = 1))
  epochs <- as.matrix(read.table(filename, nrow = 1, skip = 1))
  coal   <- read.table(filename, skip = 2)

  coal[,1] <- groups[coal[,1]+1]
  coal[,2] <- groups[coal[,2]+1]

  colnames(coal)[-c(1:2)] <- epochs
  colnames(coal)[1:2]     <- c("group1", "group2")
  coal                    <- melt(coal, id.vars = c("group1", "group2"))
  colnames(coal)[3:4]     <- c("epoch.start", "haploid.coalescence.rate")
  coal                    <- coal[order(paste(coal[,1], coal[,2], sep = "")),]

  return(coal)
}

#' Parse average rate file
#'
#' This function parses a *_avg.rate file, which stores average mutation rates inferred by Relate
#'
#' @param filename string. Filename
#' @return Returns a data frame.
#' @examples
#' read.avg_rate("example_avg.rate")
#'
read.avg_rate <- function(filename){
  rate           <- read.table(filename)
  colnames(rate) <- c("epoch.start", "mutation.rate")
  rate           <- rate[!is.na(rate[,2]),]
  return(rate)
}

#' Parse sele file
#'
#' This function parses a .sele file, which stores selection p-values inferred by Relate
#'
#' @param filename string. Filename
#' @param CHR int. Chromosome index.
#' @return Returns a data frame.
#' @examples
#' read.sele("example.sele")
#'
read.sele <- function(filename, CHR = NA){

  sele <- read.table(filename, header = T)
  colnames(sele)[-c(1:2)] <- gsub(colnames(sele)[-c(1:2)], pattern = "X", replacement = "")
  if(!is.na(CHR)) sele   <- cbind(CHR = CHR, sele)
  return(sele)

}

#' Parse lin file
#'
#' This function parses a .lin file, which stores number of lineages remaining in Relate-inferred genealogies
#'
#' @param filename string. Filename
#' @param CHR int. Chromosome index.
#' @return Returns a data frame.
#' @examples
#' read.lin("example.lin")
#'
read.lin <- function(filename, CHR = NA){

  lin                    <- read.table(filename, header = T)
  colnames(lin)[-c(1:2)] <- gsub(colnames(lin)[-c(1:2)], pattern = "X", replacement = "")
  if(!is.na(CHR)) lin    <- cbind(CHR = CHR, lin)
  return(lin)

}

#' Parse freq file
#'
#' This function parses a .freq file, which stores frequencies in Relate-inferred genealogies
#'
#' @param filename string. Filename
#' @param CHR int. Chromosome index.
#' @return Returns a data frame.
#' @examples
#' read.freq("example.freq")
#'
read.freq <- function(filename, CHR = NA){

  freq                    <- read.table(filename, header = T)
  colnames(freq)[-c(1:2)] <- gsub(colnames(freq)[-c(1:2)], pattern = "X", replacement = "")
  if(all(freq$DataFreq == 0)){
    freq$DataFreq <- NA
  }
  if(!is.na(CHR)) freq   <- cbind(CHR = CHR, freq)

  return(freq)

}

#' Parse qual file
#'
#' This function parses a .qual file, which stores summary statistics on the quality of Relate-inferred genealogies
#'
#' @param filename string. Filename
#' @param CHR int. Chromosome index.
#' @return Returns a data frame.
#' @examples
#' read.qual("example.qual")
#'
read.qual <- function(filename, CHR = NA){

  qual                   <- read.table(filename, header = T)
  if(!is.na(CHR)) qual   <- cbind(CHR = CHR, qual)
  colnames(qual)[which(colnames(qual) == "pos")] <- "BP"
  return(qual)

}


