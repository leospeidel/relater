#' Get haplotype matrix (L x N)
#'
#' This function extracts the haplotype matrix from .haps data frames
#'
#' @param haps data.table.
#' @return Returns a matrix.
#' @examples
#' haps <- read.haps(system.file("extdata/example.haps.gz", package = "relater"))
#' get_hap_matrix(haps)
#' @export

get_hap_matrix <- function(haps){

  if(ncol(haps) <= 5){
    stop("Haps file has <= 5 columns.")
  }
  hap <- as.matrix(haps[,-c(1:5)])

  if(any(is.na(hap))){
    warning("Some haplotypes coerced to NA.")
  }
  if(!all(hap == 0 | hap == 1)){
    warning("Some haplotypes are not equal to 0 or 1.")
  }
  return(hap)

}
