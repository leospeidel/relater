# polygenic selection test

#' Initialise polygenic selection test
#'
#' This function initialises the polygenic selection test.
#'
#' @param allele_ages data.table. Obtained from get.allele_ages()
#' @return Returns a list containing a data frame (quantiles) and a list of data tables (allele_ages_quantiles).
#' @examples
#' PolyTest_Init(allele_ages)
#' @export

PolyTest_Init <- function(allele_ages){
  allele_ages <- allele_ages[!is.na(allele_ages$pvalue),]
  allele_ages <- allele_ages[allele_ages$pvalue < 0,]

  ##divide SNPs by frequency quantile
  quantiles                    <- stats::quantile(allele_ages$DAF, probs=seq(0,1,0.05))
  quantiles                    <- unique(quantiles)
  quantiles[length(quantiles)] <- max(allele_ages$DAF+1)
  freq_quantiles               <- data.table(ID = allele_ages$ID, DAF = cut(allele_ages$DAF,  quantiles, include.lowest=TRUE))
  allele_ages_quantiles        <- data.table::split(data.table(ID = allele_ages$ID, pvalue = allele_ages$pvalue, DAF = allele_ages$DAF), with(freq_quantiles, DAF))
  quantiles                    <- quantiles[-1]

  return(list(quantiles = quantiles, allele_ages_quantiles = allele_ages_quantiles))
}

#' Resample selection p-values
#'
#' Given DAF, resample 20 selection p-values of SNPs selected uniformly at random with matching DAF
#'
#' @param DAF int. Derived allele frequency
#' @param quantiles data frame. Obtained from PolyTest_Init.
#' @param allele_ages_quantiles List of data tables. Obtained from PolyTest_Init.
#' @return Returns a numeric 1x20 matrix.
#' @examples
#' PolyTest_ResampleSNPs(DAF, quantiles, allele_ages_quantiles)
#' @export
#'
PolyTest_ResampleSNPs <- function(DAF, quantiles, allele_ages_quantiles){

  sample      <- data.frame(DAF = DAF, n = 1)

  #resample with same frequencies as in sample
  resampled <- mapply(function(i){
    #get the right quantile and then match exactly within the right quantile
    foo <- allele_ages_quantiles[[min(which(sample$DAF[i] <= quantiles))]][DAF == as.numeric(as.character(sample$DAF[i])), "pvalue"][sample(.N, sample$n[i] * 20,replace = T)]
    return(foo)
  },1:dim(sample)[1])

  resampled <- unlist(resampled)
  return(t(as.numeric(as.matrix(resampled))))

}

#' Test for evidence of polygenic selection
#'
#' Tests for evidence of polygenic selection by matching SNPs by DAF and using a one-sided Wilcoxon rank sum test.
#'
#' @param DAFs 1d array. Derived allele frequencies of trait associations
#' @param pvalues 1d array. Selection pvalues of trait associations
#' @param quantiles data frame. Obtained from PolyTest_Init.
#' @param allele_ages_quantiles List of data tables. Obtained from PolyTest_Init.
#' @return Returns a pvalue.
#' @examples
#' PolyTest(DAFs, pvalues, quantiles, allele_ages_quantiles)
#'
#' # Example analysis:
#' # read files
#' mut         <- read.mut("./inst/extdata/example.mut.gz")
#' sele        <- read.sele("./inst/extdata/example.sele.gz")
#' freq        <- read.freq("./inst/extdata/example.freq.gz")
#' qual        <- read.qual("./inst/extdata/example.qual.gz")
#'
#' # Obtain allele_ages data table
#' allele_ages <- get.allele_ages(mut, freq, sele)
#' allele_ages <- filter.allele_ages(allele_ages, qual)
#'
#' ######### Polygenic selection test #########
#'
#' # Initialise
#' quant       <- PolyTest_Init(allele_ages)
#'
#' # Make a fake polygenic trait
#' df          <- allele_ages[!is.na(allele_ages$pvalue),]
#' df          <- df[sample(1:nrow(df), 50, replace = F),]
#'
#' # Run polygenic test
#' PolyTest(df$DAF, df$pvalue, quant$quantiles, quant$allele_ages_quantiles)
#' @export

PolyTest <- function(DAFs, pvalues, quantiles, allele_ages_quantiles){

  if(length(DAFs) != length(pvalues)){
    stop("DAFs and pvalues have different length.")
  }

  p <- mapply(function(k){

    resampled <- numeric(0)
    #for each SNP, resample 20 pvalues with matched DAF
    for(i in 1:length(DAFs)){
      resampled <- c(resampled, PolyTest_ResampleSNPs(DAFs[i], quantiles, allele_ages_quantiles))
    }

    #use wilcox test to test if pvalues are shifted compared to random SNPs
    stats::wilcox.test(pvalues, resampled, alternative = "less", mu = 0)$p.value

  },k = 1:20)

  return(mean(p))

}
