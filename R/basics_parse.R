#' Parse haps file
#'
#' This function uses fread to parse a .haps file, which is part of the haps/sample output file format of SHAPEIT2.
#'
#' @param filename string. Filename
#' @param ... Any other parameters for fread
#' @return Returns a data table.
#' @examples
#' read.haps(system.file("extdata/example.haps.gz", package = "relater"))
#' @export

read.haps <- function(filename, ...){
  haps <- data.table::fread(filename, ...)
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
#' read.sample(system.file("extdata/example.sample.gz", package = "relater"))
#' @export

read.sample <- function(filename, ...){
  sam <- data.table::fread(filename, ...)
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
#' read.mut(system.file("extdata/example.mut.gz", package = "relater"))
#' @export
#'
read.mut <- function(filename, CHR = NA, ...){
  mut <- data.table::fread(filename, ...)
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
#' read.coal(system.file("extdata/example.coal.gz", package = "relater"))
#' @export
#'
read.coal <- function(filename){

  groups <- as.matrix(utils::read.table(filename, nrow = 1))
  epochs <- as.matrix(utils::read.table(filename, nrow = 1, skip = 1))
  coal   <- utils::read.table(filename, skip = 2)

  coal[,1] <- groups[coal[,1]+1]
  coal[,2] <- groups[coal[,2]+1]

  colnames(coal)[-c(1:2)] <- epochs
  colnames(coal)[1:2]     <- c("group1", "group2")
  coal                    <- reshape2::melt(coal, id.vars = c("group1", "group2"))
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
#' read.avg_rate(system.file("extdata/example_avg.rate.gz", package = "relater"))
#' @export
#'
read.avg_rate <- function(filename){
  rate           <- utils::read.table(filename)
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
#' read.sele(system.file("extdata/example.sele.gz", package = "relater"))
#' @export
#'
read.sele <- function(filename, CHR = NA){

  sele <- utils::read.table(filename, header = T)
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
#' read.lin(system.file("extdata/example.lin.gz", package = "relater"))
#' @export
read.lin <- function(filename, CHR = NA){

  lin                    <- utils::read.table(filename, header = T)
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
#' read.freq(system.file("extdata/example.freq.gz", package = "relater"))
#' @export
read.freq <- function(filename, CHR = NA){

  freq                    <- utils::read.table(filename, header = T)
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
#' read.qual(system.file("extdata/example.qual.gz", package = "relater"))
#' @export
read.qual <- function(filename, CHR = NA){

  qual                   <- utils::read.table(filename, header = T)
  if(!is.na(CHR)) qual   <- cbind(CHR = CHR, qual)
  colnames(qual)[which(colnames(qual) == "pos")] <- "BP"
  return(qual)

}


